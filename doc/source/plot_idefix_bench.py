import matplotlib.pyplot as plt
import json

def do_plot(title, bench_file, gpumodel):
    with open(bench_file, 'r') as f:
        benches = json.load(f)

    select = [bench for bench in benches if bench['gpumodel'] == gpumodel][-1]
    xs = [r['nbgpu'] for r in select['results']]
    ys = [r['cell_updates'] for r in select['results']]

    plt.xscale("log", base=2)
    plt.plot(xs, ys)
    plt.ylim(0,max(ys)*1.1)
    plt.xlabel("Number of GPUs")
    plt.ylabel("Performance (cells / second / GPU)")
    plt.title(title)
