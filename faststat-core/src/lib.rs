use std::cmp::Ordering as CmpOrdering;
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(target_os = "windows")]
fn raw_nanotime() -> u64 {
    use windows_sys::Win32::Foundation::FILETIME;
    use windows_sys::Win32::System::SystemInformation::GetSystemTimeAsFileTime;
    const DELTA_EPOCH_IN_SECS: u64 = 11_644_473_600;
    unsafe {
        let mut ft = FILETIME {
            dwLowDateTime: 0,
            dwHighDateTime: 0,
        };
        GetSystemTimeAsFileTime(&mut ft);
        let ticks = ((ft.dwHighDateTime as u64) << 32) | (ft.dwLowDateTime as u64);
        (ticks - DELTA_EPOCH_IN_SECS * 10_000_000) * 100
    }
}

#[cfg(all(unix, not(target_os = "macos")))]
fn raw_nanotime() -> u64 {
    use libc::{clock_gettime, timespec, CLOCK_REALTIME};
    unsafe {
        let mut ts: timespec = std::mem::zeroed();
        if clock_gettime(CLOCK_REALTIME, &mut ts) == 0 {
            (ts.tv_sec as u64) * 1_000_000_000 + ts.tv_nsec as u64
        } else {
            0
        }
    }
}

#[cfg(target_os = "macos")]
fn raw_nanotime() -> u64 {
    use libc::{gettimeofday, timeval};
    unsafe {
        let mut tv: timeval = std::mem::zeroed();
        if gettimeofday(&mut tv, std::ptr::null_mut()) == 0 {
            (tv.tv_sec as u64) * 1_000_000_000 + (tv.tv_usec as u64) * 1_000
        } else {
            0
        }
    }
}

static NANOTIME_OVERRIDE: AtomicU64 = AtomicU64::new(0);

#[derive(Clone)]
struct P2Percentile {
    percentile: u16,
    val: f64,
    n: f64,
}

#[derive(Clone)]
struct Bucket {
    max: f64,
    count: u64,
}

#[derive(Clone, Default)]
struct DataPoint {
    value: f64,
    nanotime: u64,
}

#[derive(Clone)]
struct ExpoAvg {
    val: f64,
    alpha: f64,
}

#[derive(Clone)]
struct WindowCount {
    num_windows: usize,
    window_size: u64,
    counts: Vec<u32>,
}

#[derive(Clone)]
pub struct StatsLite {
    n: u64,
    mean: f64,
    m2: f64,
    m3: f64,
    m4: f64,
    min: f64,
    max: f64,
    sum_of_logs: f64,
    sum_of_inv: f64,
}

impl StatsLite {
    pub fn new() -> Self {
        StatsLite {
            n: 0,
            mean: 0.0,
            m2: 0.0,
            m3: 0.0,
            m4: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sum_of_logs: 0.0,
            sum_of_inv: 0.0,
        }
    }

    fn update_moments(&mut self, x: f64) {
        let n = self.n as f64;
        let delta = x - self.mean;
        let delta_n = delta / n;
        let delta_m2 = delta * delta_n * (n - 1.0);
        let delta_m3 = delta_m2 * delta_n * (n - 2.0);
        let delta_m4 = delta_m2 * delta_n * delta_n * (n * (n - 3.0) + 3.0);
        self.mean += delta_n;
        self.m4 += delta_m4 + delta_n * (6.0 * delta_n * self.m2 - 4.0 * self.m3);
        self.m3 += delta_m3 + delta_n * 3.0 * self.m2;
        self.m2 += delta_m2;
    }

    pub(crate) fn add_internal(&mut self, x: f64) {
        self.n += 1;
        if self.n == 1 {
            self.min = x;
            self.max = x;
        }
        if x <= self.min {
            self.min = x;
        }
        if x >= self.max {
            self.max = x;
        }
        self.sum_of_logs += if x > 0.0 { x.ln() } else { f64::NAN };
        self.sum_of_inv += if x > 0.0 { 1.0 / x } else { f64::NAN };
    }

    pub fn add(&mut self, x: f64) {
        self.add_internal(x);
    }

    pub fn variance(&self) -> Option<f64> {
        if self.n < 2 { None } else { Some(self.m2 / (self.n as f64 - 1.0)) }
    }

    pub fn n(&self) -> u64 { self.n }
    pub fn mean(&self) -> f64 { self.mean }
    pub fn m2(&self) -> f64 { self.m2 }
    pub fn m3(&self) -> f64 { self.m3 }
    pub fn m4(&self) -> f64 { self.m4 }
    pub fn min(&self) -> f64 { self.min }
    pub fn max(&self) -> f64 { self.max }
    pub fn sum_of_logs(&self) -> f64 { self.sum_of_logs }
    pub fn sum_of_inv(&self) -> f64 { self.sum_of_inv }
}

#[derive(Clone)]
pub struct Stats {
    core: StatsLite,
    mintime: u64,
    maxtime: u64,
    lasttime: u64,
    percentiles: Vec<P2Percentile>,
    buckets: Vec<Bucket>,
    expo_avgs: Vec<ExpoAvg>,
    window_avg: f64,
    lastn: Vec<DataPoint>,
    topn: Vec<DataPoint>,
    window_counts: Vec<WindowCount>,
    interval: Option<Box<Stats>>,
}

impl Stats {
    fn nanotime() -> u64 {
        let val = NANOTIME_OVERRIDE.load(Ordering::Relaxed);
        if val != 0 { val } else { raw_nanotime() }
    }

    fn update_moments(&mut self, x: f64) {
        self.core.update_moments(x);
    }

    fn insert_percentile_sorted(&mut self, mut x: f64) {
        let num = self.core.n.min(self.percentiles.len() as u64) as usize;
        if num == 0 { return; }
        for i in 0..num - 1 {
            if x < self.percentiles[i].val {
                std::mem::swap(&mut x, &mut self.percentiles[i].val);
            }
        }
        self.percentiles[num - 1].val = x;
    }

    fn p2_update_point(l_v: f64, l_n: f64, cur: &mut P2Percentile, r_v: f64, r_n: f64, n: f64) {
        let percentile = cur.percentile as f64 / 65536.0;
        let diff = (n - 1.0) * percentile + 1.0 - cur.n;
        let d = if diff >= 1.0 { 1.0 } else if diff <= -1.0 { -1.0 } else { return };
        let c_v = cur.val;
        let new_val = if l_n < cur.n + d && cur.n + d < r_n {
            let val = c_v
                + (d / (r_n - l_n))
                    * ((cur.n - l_n + d) * (r_v - c_v) / (r_n - cur.n)
                        + (r_n - cur.n - d) * (c_v - l_v) / (cur.n - l_n));
            if l_v >= val || r_v <= val {
                if d == 1.0 {
                    c_v + (r_v - c_v) / (r_n - cur.n)
                } else {
                    c_v - (l_v - c_v) / (cur.n - l_n)
                }
            } else {
                val
            }
        } else if d == 1.0 {
            c_v + (r_v - c_v) / (r_n - cur.n)
        } else {
            c_v - (l_v - c_v) / (cur.n - l_n)
        };
        cur.val = new_val;
        cur.n += d;
    }

    fn update_percentiles(&mut self, x: f64) {
        if self.percentiles.is_empty() { return; }
        if self.core.n <= self.percentiles.len() as u64 {
            self.insert_percentile_sorted(x);
            return;
        }
        let n = self.core.n as f64;
        let last = self.percentiles.len() - 1;
        if x < self.percentiles[0].val {
            self.percentiles[0].val = x;
        } else if x >= self.percentiles[last].val {
            self.percentiles[last].val = x;
            self.percentiles[last].n = n;
        } else {
            for i in 1..=last {
                if x < self.percentiles[i].val {
                    for j in i..=last {
                        self.percentiles[j].n += 1.0;
                    }
                    break;
                }
            }
        }
        let mut prev_v = self.core.min;
        let mut prev_n = 0.0;
        for i in 0..last {
            let (r_v, r_n) = {
                let nxt = &self.percentiles[i + 1];
                (nxt.val, nxt.n)
            };
            let cur = &mut self.percentiles[i];
            Self::p2_update_point(prev_v, prev_n, cur, r_v, r_n, n);
            prev_v = cur.val;
            prev_n = cur.n;
        }
        let cur = &mut self.percentiles[last];
        Self::p2_update_point(prev_v, prev_n, cur, self.core.max, n, n);
    }

    fn update_buckets(&mut self, x: f64) {
        for b in &mut self.buckets {
            if x < b.max {
                b.count += 1;
                break;
            }
        }
    }

    fn update_expo_avgs(&mut self, x: f64) {
        for ea in &mut self.expo_avgs {
            ea.val = x * ea.alpha + ea.val * (1.0 - ea.alpha);
        }
    }

    fn update_lastn(&mut self, x: f64) {
        if self.lastn.is_empty() { return; }
        let len = self.lastn.len();
        let idx = ((self.core.n() - 1) as usize) & (len - 1);
        self.window_avg -= self.lastn[idx].value / len as f64;
        self.window_avg += x / len as f64;
        self.lastn[idx].value = x;
        self.lastn[idx].nanotime = self.lasttime;
    }

    fn rezero_window_counts(&mut self, t: u64) {
        for wc in &mut self.window_counts {
            let last_window = self.lasttime / wc.window_size;
            let cur_window = t / wc.window_size;
            if last_window == cur_window { continue; }
            let diff = cur_window.saturating_sub(last_window);
            if diff as usize >= wc.num_windows {
                wc.counts.fill(0);
            } else {
                for j in 1..=diff {
                    let idx = ((last_window + j) & (wc.num_windows as u64 - 1)) as usize;
                    wc.counts[idx] = 0;
                }
            }
        }
    }

    fn update_window_counts(&mut self, t: u64) {
        self.rezero_window_counts(t);
        for wc in &mut self.window_counts {
            let idx = ((t / wc.window_size) & (wc.num_windows as u64 - 1)) as usize;
            wc.counts[idx] += 1;
        }
    }

    fn update_topn(&mut self, x: f64, t: u64) {
        if self.topn.is_empty() { return; }
        self.topn.push(DataPoint { value: x, nanotime: t });
        self.topn.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap_or(CmpOrdering::Equal));
        if self.topn.len() > self.num_top() { self.topn.pop(); }
    }

    fn num_top(&self) -> usize { self.topn.capacity() }
    fn num_prev(&self) -> usize { self.lastn.capacity() }
}

impl Stats {
    pub fn with_options(
        buckets: &[f64],
        num_prev: usize,
        percentiles: &[f64],
        interval: Option<Box<Stats>>,
        expo_avgs: &[f64],
        window_counts: &[(usize, u64)],
        num_top: usize,
    ) -> Self {
        let percentiles = percentiles
            .iter()
            .enumerate()
            .map(|(i, p)| P2Percentile {
                percentile: ((*p * 65536.0).round() as u16),
                val: 0.0,
                n: (i + 1) as f64,
            })
            .collect();
        let buckets = buckets
            .iter()
            .map(|b| Bucket { max: *b, count: 0 })
            .collect();
        let expo_avgs = expo_avgs
            .iter()
            .map(|a| ExpoAvg { val: 0.0, alpha: *a })
            .collect();
        let window_counts = window_counts
            .iter()
            .map(|(n, size)| WindowCount {
                num_windows: *n,
                window_size: *size,
                counts: vec![0; *n],
            })
            .collect();
        let lastn = vec![DataPoint::default(); num_prev];
        let topn: Vec<DataPoint> = Vec::with_capacity(num_top);
        Stats {
            core: StatsLite::new(),
            mintime: 0,
            maxtime: 0,
            lasttime: 0,
            percentiles,
            buckets,
            expo_avgs,
            window_avg: 0.0,
            lastn,
            topn,
            window_counts,
            interval,
        }
    }

    pub fn new() -> Self {
        Self::with_options(&[], 0, &[], None, &[], &[], 0)
    }

    fn _add(&mut self, x: f64, t: u64) {
        self.lasttime = t;
        self.core.add_internal(x);
        if self.core.n == 1 {
            self.mintime = t;
            self.maxtime = t;
        }
        if x <= self.core.min {
            self.mintime = t;
        }
        if x >= self.core.max {
            self.maxtime = t;
        }
        self.update_percentiles(x);
        self.update_buckets(x);
        self.update_expo_avgs(x);
        self.update_lastn(x);
        self.update_window_counts(t);
        self.update_topn(x, t);
    }

    pub fn add(&mut self, x: f64) {
        let t = Self::nanotime();
        if let Some(ref mut interval) = self.interval {
            if self.lasttime != 0 {
                let diff = (t - self.lasttime).max(1);
                interval._add(diff as f64, t);
            }
        }
        self._add(x, t);
    }

    pub fn end(&mut self, start: u64) {
        let end = Self::nanotime();
        if let Some(ref mut interval) = self.interval {
            if self.lasttime != 0 {
                let diff = (end - self.lasttime).max(1);
                interval._add(diff as f64, end);
            }
        }
        self._add((end - start) as f64, end);
    }

    pub fn tick(&mut self) {
        let t = Self::nanotime();
        if self.lasttime != 0 {
            let diff = (t - self.lasttime) as f64;
            self._add(diff, t);
        } else {
            self.lasttime = t;
        }
    }

    pub fn get_percentiles(&self) -> Vec<(f64, f64)> {
        self.percentiles
            .iter()
            .map(|p| {
                let key = ((p.percentile as f64) / 65536.0 * 10000.0).floor() / 10000.0;
                (key, p.val)
            })
            .collect()
    }

    pub fn get_buckets(&self) -> Vec<(f64, u64)> {
        let mut leftover = self.core.n();
        let mut out = Vec::new();
        for b in &self.buckets {
            leftover -= b.count;
            out.push((b.max, b.count));
        }
        out.push((f64::NAN, leftover));
        out
    }

    pub fn get_expo_avgs(&self) -> Vec<(f64, f64)> {
        self.expo_avgs.iter().map(|ea| (ea.alpha, ea.val)).collect()
    }

    pub fn get_prev(&self, offset: usize) -> Option<(u64, f64)> {
        if self.lastn.is_empty() { return None; }
        let len = self.num_prev();
        let idx = (((self.core.n() - 1) as usize + (len - offset)) & (len - 1)) as usize;
        let dp = &self.lastn[idx];
        Some((dp.nanotime, dp.value))
    }

    pub fn get_top_n(&self) -> Vec<(f64, u64)> {
        self.topn.iter().map(|dp| (dp.value, dp.nanotime)).collect()
    }

    pub fn get_window_counts(&self) -> Vec<(u64, Vec<u32>)> {
        let t = Self::nanotime();
        let wc_clone = self.window_counts.clone();
        let mut temp = Stats { window_counts: wc_clone.clone(), ..self.clone() };
        temp.rezero_window_counts(t);
        temp.window_counts
            .into_iter()
            .map(|wc| {
                let cur_window = t / wc.window_size;
                let mut items = Vec::with_capacity(wc.num_windows);
                for j in 0..wc.num_windows {
                    let idx = ((cur_window - j as u64) & (wc.num_windows as u64 - 1)) as usize;
                    items.push(wc.counts[idx]);
                }
                (wc.window_size, items)
            })
            .collect()
    }

    pub fn n(&self) -> u64 { self.core.n() }
    pub fn mean(&self) -> f64 { self.core.mean() }
    pub fn m2(&self) -> f64 { self.core.m2() }
    pub fn m3(&self) -> f64 { self.core.m3() }
    pub fn m4(&self) -> f64 { self.core.m4() }
    pub fn variance(&self) -> Option<f64> { self.core.variance() }
    pub fn min(&self) -> f64 { self.core.min() }
    pub fn max(&self) -> f64 { self.core.max() }
    pub fn sum_of_logs(&self) -> f64 { self.core.sum_of_logs() }
    pub fn sum_of_inv(&self) -> f64 { self.core.sum_of_inv() }
    pub fn lasttime(&self) -> u64 { self.lasttime }
}

pub fn nanotime() -> u64 { Stats::nanotime() }

pub fn nanotime_override(t: u64) { NANOTIME_OVERRIDE.store(t, Ordering::Relaxed); }

