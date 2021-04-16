f=open("degeneslist.txt","r")

L=[]
Lines = f.readlines()
count = 0
# Strips the newline character
for line in Lines:
    wor=line
    if("(" in line):
        # wor= line[line.index("(")+1:line.index(")")]
        wor=line[ line.rfind("(")+1 : line.rfind(")")]
        count+=1
        L.append(wor+"\n")

print(L)
print(len(L))
file1 = open('out2.txt', 'w')
file1.writelines(L)
file1.close()