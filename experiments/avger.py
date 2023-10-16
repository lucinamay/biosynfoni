from sys import argv
lines  =[]
with open(argv[1], 'r') as f:
    for line in f:
        lines.append(float(line.strip()))

print(sum(lines)/len(lines))

