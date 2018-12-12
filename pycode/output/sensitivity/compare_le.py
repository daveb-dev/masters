'''
This program was originally used to determine the difference in times when running two different 
versions of the forward model with linear elasticity. The average difference was calculated. 
It can now be used to compare any two files of the same format with running times every nth line.
'''

if(len(sys.argv) != 4):
    print("Wrong number of inputs. USAGE: python compare_le.py file1 file2 n")
    quit()
file1 = sys.argv[1]
file2 = sys.argv[2]
n     = int(sys.argv[3])

le1 = open(file1)
le2 = open(file2)
lines1 = le1.readlines()
lines2 = le2.readlines()
i = 0
diff=[]
for l in lines1:
    if ((i+1)%n) == 0:
        diff.append(float(lines1[i])-float(lines2[i]))
        print('diff = '+str(diff[-1]))
    i += 1
print('average diff = '+str(sum(diff)/float(len(diff))))
