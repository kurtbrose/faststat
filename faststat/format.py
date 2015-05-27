import math
import json
 

def stat2json(stat):
    prev = [stat.get_prev(i) for i in range(stat.num_prev)]
    # timestamp = 0 means not yet valid
    prev = [(p[0] / 1e6, p[1]) for p in prev if p[0]]
    return json.dumps({
        "n": stat.n,
        "mean": stat.mean,
        "max": stat.max,
        "min": stat.min,
        "percentiles": stat.get_percentiles(),
        "prev": prev
    })


def stat2html(stat):
    return TEMPLATE.replace('"==THE_STAT=="', stat2json(stat))


def si_round(val):
    '''
    round to a "scientific notation" tuple of (factor, exponent)
    such that 1 < factor < 1000, and factor * 10 ** exponent == val
    '''
    if val < 0:
        neg = True
        val = -val
    elif val == 0:
        return 0, 0
    else:
        neg = False
    exp = math.log(val) / math.log(1000)
    if exp < 0:
        exp = int(exp) - 1
    else:
        exp = int(exp)
    val = val / 1000.0 ** exp
    if neg:
        val = -val
    return val, 3 * exp

 
def si_format(val):
    val, exp = si_round(val)
    if exp:
        exps = _SCALES.get(exp) or 'E' + str(exp)
    else:
        exps = ''
    if val >= 100 or val <= -100:
        return '{0:0.0f}{1}'.format(val, exps)
    if val >= 10 or val <= -10:
        return '{0:0.1f}{1}'.format(val, exps)
    return '{0:0.2f}{1}'.format(val, exps)


def sib_round(val):
    '''
    round to a binary SI tuple of (factor, exponent)
    such that 1 < factor < 1024, and factor * 1024 ** exponent == val
    '''
    if val < 0:
        neg = True
        val = -val
    elif val == 0:
        return 0, 0
    else:
        neg = False
    exp = math.log(val) / math.log(1024)
    if exp < 0:
        exp = int(exp) - 1
    else:
        exp = int(exp)
    val = val / 1024.0 ** exp
    if neg:
        val = -val
    return val, exp


def sib_format(val):
    val, exp = sib_round(val)
    if exp < 0 or exp > len(_BSCALES):
        raise ValueError("{0} out of format range", val)
    exps = _BSCALES[exp]
    if val >= 100 or val <= -100:
        return '{0:0.0f}{1}'.format(val, exps)
    if val >= 10 or val <= -10:
        return '{0:0.1f}{1}'.format(val, exps)
    return '{0:0.2f}{1}'.format(val, exps)


_SCALES = {
    # Peta, Tera, Giga, Mega, kilo
    15: 'P', 12: 'T', 9: 'G', 6: 'M', 3: 'k',
    -3: 'm', -6: 'u', -9: 'n', -12: 'p', -15: 'f'
    # milli, micro, nano, pico, femto
}


_BSCALES = ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi']


JAVASCRIPT_HTML_HEAD = '''
  <script src="http://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js"></script>
  <script src="http://dimplejs.org/dist/dimple.v2.1.2.min.js"></script>
  <script type="text/javascript">
  // histogram is computed from percentiles by taking the region over which the percentile
  // applies, dividing by the fraction of points that occur in that region
  /// e.g. if 25% = 0.18, then X = min, Y = 0.25 / (0.18 - min)
    function faststat_histogram_chart(faststat, container) {
      var svg = dimple.newSvg(container, 590, 400);
      // var avg_density = faststat.n / (faststat.max - faststat.min);
      var width = faststat.max - faststat.min
      var chart = new dimple.chart(svg);
      var x = chart.addMeasureAxis("x", "Value");
      chart.addMeasureAxis("y", "Density");
      var histogram = chart.addSeries("percentile", dimple.plot.line);
      histogram.interpolation = 'step-before';
      histogram.data = [
        //the value is the minimum, the density is percentage divided by width
        {"percentile": "min", "Value": faststat.min, "Density": 0},
        {"percentile": 0, "Value": faststat.min, "Density": (0.25 - 0) / (faststat.percentiles[0.25] - faststat.min)},
        {"percentile": 25, "Value": faststat.percentiles[0.25], "Density": (0.5 - 0.25) * width / (faststat.percentiles[0.5] - faststat.percentiles[0.25])},
        {"percentile": 50, "Value": faststat.percentiles[0.5], "Density": (0.75 - 0.5) * width / (faststat.percentiles[0.75] - faststat.percentiles[0.5])},
        {"percentile": 75, "Value": faststat.percentiles[0.75], "Density": (0.9 - 0.75) * width / (faststat.percentiles[0.9] - faststat.percentiles[0.75])},
        {"percentile": 90, "Value": faststat.percentiles[0.9], "Density": (0.95 - 0.9) * width / (faststat.percentiles[0.95] - faststat.percentiles[0.9])},
        {"percentile": 95, "Value": faststat.percentiles[0.95], "Density": (0.99 - 0.95) * width / (faststat.percentiles[0.99] - faststat.percentiles[0.95])},
        {"percentile": 99, "Value": faststat.percentiles[0.99], "Density": (1.0 - 0.99) * width / (faststat.max - faststat.percentiles[0.99])},
        {"percentile": 100, "Value": faststat.max, "Density": 0}
      ];
      for(var i=0; i<histogram.data.length; i++) {
        if(histogram.data[i].Density === NaN) {
          histogram.data[i].Density = 1000;
        }
      }
      var vertical = chart.addPctAxis("y", "Vertical");
      vertical.hidden = true;
      var mean = chart.addSeries("mean", dimple.plot.area, [x, vertical]);
      mean.data = [
        {"Value": faststat.mean, "Vertical": 0, "mean": 0},
        {"Value": faststat.mean, "Vertical": 1, "mean": 1},
      ];
      chart.addLegend(60, 10, 500, 20, "right");
      chart.draw();
    }

    function faststat_time_chart(faststat, container) {
      var svg = dimple.newSvg(container, 590, 400);
      var chart = new dimple.chart(svg);
      var x = chart.addMeasureAxis("x", "timestamp");
      var y = chart.addMeasureAxis("y", "Value");
      var recent = chart.addSeries("index", dimple.plot.line);
      var recent_data = [];
      for(var i=0; i<faststat.prev.length; i++) {
        recent_data.push({"index": i, "timestamp": faststat.prev[i][0] - faststat.prev[0][0], "Value": faststat.prev[i][1]});
      }
      console.log(recent_data);
      recent.data = recent_data;
      chart.draw();
    }
  </script>
'''

TEMPLATE = '''
<html>
<body>
<div id="histogram_container">
<div id="time_container">
  <!==JAVASCRIPT==>
  <script type="text/javascript">
    var THE_STAT = "==THE_STAT==";
    faststat_histogram_chart(THE_STAT, "#histogram_container");
    faststat_time_chart(THE_STAT, "#time_container");
  </script>
</div>
</body>
</html>'''.replace('<!==JAVASCRIPT==>', JAVASCRIPT_HTML_HEAD)
