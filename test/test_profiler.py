import faststat.profiler as profiler


def busy_wait():
    for _ in range(1000):
        pass


def test_function_stats():
    p = profiler.Profiler(func_stats=True, line_stats=False)
    p.start()
    busy_wait()
    p.stop()
    key = (
        busy_wait.__code__.co_filename,
        busy_wait.__name__,
        busy_wait.__code__.co_firstlineno,
    )
    assert key in p.func_data
    stat = p.func_data[key]
    assert stat.n >= 1
    assert stat.mean > 0


def test_line_stats():
    p = profiler.Profiler(func_stats=False, line_stats=True)
    p.start()
    busy_wait()
    p.stop()
    file = busy_wait.__code__.co_filename
    assert any(k[0] == file for k in p.line_data)

