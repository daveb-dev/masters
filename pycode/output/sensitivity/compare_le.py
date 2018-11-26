le1 = open('le_log.txt')
le2 = open('le_log2.txt')
lines1 = le1.readlines()
lines2 = le2.readlines()
i = 0
diff=[]
for l in lines1:
    if ((i+1)%3) == 0:
        diff.append(float(lines1[i])-float(lines2[i]))
        print('diff = '+str(diff[-1]))
    i += 1
print('average diff = '+str(sum(diff)/float(len(diff))))
