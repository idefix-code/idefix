import matplotlib.pyplot as plt
import json

def do_plot(title, bench_file, gpumodel):
    with open(bench_file, 'r') as f:
        benches = json.load(f)

    select = [bench for bench in benches if bench['gpumodel'] == gpumodel][-1]
    xs = [r['nbgpu'] for r in select['results']]
    ys = [r['cell_updates'] for r in select['results']]
    
    plt.xscale("log")
    plt.plot(xs, ys)
    plt.title(title)
