import datetime
file = open("./cleanestcalls.csv", "r")
lines = file.readlines()
orig = datetime.datetime.strptime("03/03/20", "%m/%d/%y")
for x,i in enumerate(lines):
    args = i.split(" ")

    exp = datetime.datetime.strptime(args[1], "%m/%d/%y")
    args[1] = str(round((exp-orig).days / 365, 5))
    lines[x] = " ".join(args)
file.close()
file = open("./cleanestcalls.csv", "w")
file.write("".join(lines))
file.close()
