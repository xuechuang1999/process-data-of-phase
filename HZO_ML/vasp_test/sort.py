from glob import glob
lis=glob('./POSCAR-*')
def dis(coord1,coord2):
    a=0
    for i in range(3):
        a+=min(coord1[i]-coord2[i],abs((coord1[i]-coord2[i])%1-1))**2
    return a**0.5
for k in lis:
    with open(k,'r') as file:
        content=file.readlines()
        hf=content[8:12]
        hf1=sorted(hf,key=lambda x:dis([float(i) for i in x.split()],[0.0,0.0,0.0]))[0]
        hf2=sorted(hf,key=lambda x:dis([float(i) for i in x.split()],[0.0,0.5,0.5]))[0]
        hf3=sorted(hf,key=lambda x:dis([float(i) for i in x.split()],[0.5,0.0,0.5]))[0]
        hf4=sorted(hf,key=lambda x:dis([float(i) for i in x.split()],[0.5,0.5,0.0]))[0]
        o=sorted(content[12:20],key=lambda x:[float(i)<=0.5 for i in x.split()])
        newcontent=content[:8]+[hf1,hf2,hf3,hf4]+o
        newfile=open(k+'.sorted','w')
        newfile.writelines(newcontent)
        newfile.close()
