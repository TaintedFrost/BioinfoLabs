text="actgaaactgttggcatgct"
f={}

for char in text:
    if char.isalpha():
        f[char]=f.get(char, 0) + 1

for char, count in f.items():
    print(char, ":", count)

def find_f(text, length):
    f={}
    for i in range(len(text) - length + 1):
        combination=text[i:i + length]

        if combination.isalpha():
            f[combination]=f.get(combination, 0) + 1

    return f

two_letters=find_f(text, 2)
print("\n2 letter")
for combination, count in sorted(two_letters.items(),key=lambda x: x[1],reverse=True):
    print(combination, ":", count)


three_letters=find_f(text, 3)

print("\n3 letter")
for combination, count in sorted(three_letters.items(),key=lambda x: x[1],reverse=True):
    print(combination, ":", count)