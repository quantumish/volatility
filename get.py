import yahoo_fin.options as ops
import statistics
df = ops.get_calls("ual")
print(len(df["Implied Volatility"]), "calls")
print("Mean σ:",  statistics.mean([float(i[:-1]) for i in df["Implied Volatility"]])/100)
