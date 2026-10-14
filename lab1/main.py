text="actgaaactgttggcatgct"
f={}

for char in text:
    if char.isalpha():
        f[char]=f.get(char,0)+1

for char, count in f.items():
    print(char,":",count)