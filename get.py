import yahoo_fin.options as ops
import statistics
df = ops.get_calls("ual")
print(statistics.mean([float(i[:-1]) for i in df["Implied Volatility"]]))
