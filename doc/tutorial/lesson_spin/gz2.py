import os,string,sys,math
#First find the coordinates of atoms
f=open('first_round','w')
f.write("tspin_2o_DEN\n")
f.write("1\n-1\n13\n")
f.close()
start_atomic=0
atomic_positions=[]
for line in os.popen("../../cut3d< first_round").readlines():
        if(start_atomic==1):
                if(len(line)<3):
                        start_atomic=0
                        break
                coord=string.split(line)[1:]
                atomic_positions.append(coord)
        if(string.find(line,"Atomic positions")>0):
                start_atomic=1
print "Number of atoms = ", len(atomic_positions)
print "Atomic coordinates"
print atomic_positions
num_atoms=len(atomic_positions)
cube_side=3.3
cube_step=0.1
radius=cube_side/2.0
npts=int(cube_side/cube_step)
print "number of integration points:",npts
g=open('data','w')
integral=[]
for iatom in range(num_atoms):
        npts_integral=0
        sum=0
        f=open('sphere','w')
        f.write("tspin_2o_DEN\n")
        f.write("1\n-1\n1\n1\n")
        f.write("0.0 0.0 .0 \n")
        x0=string.atof(atomic_positions[iatom][0])-cube_side/2.0
        y0=string.atof(atomic_positions[iatom][1])-cube_side/2.0
        z0=string.atof(atomic_positions[iatom][2])-cube_side/2.0
        print "Treating  atom #",iatom
        for i in range(npts):
                for j in range(npts):
                        for k in range(npts):
                                #print i,j,k
                                x=x0+i*cube_step
                                y=y0+j*cube_step
                                z=z0+k*cube_step
                                #if((x-x0)**2+(y-y0)**2+(z-z0)**2<radius**2):
                                f.write('1'+"\n"+'1'+'\n'+'1'+'\n')
                                f.write(`x`+' '+`y`+' '+`z`+'\n')
        f.write('2\n')
        f.close()
        #sys.exit()
        for line in os.popen("../../cut3d< sphere").readlines():
                if(string.find(line,"Calculated")>0):
                        npts_integral=npts_integral+1
                        g.write(string.split(line)[2]+'\n')
                        sum=sum+string.atof(string.split(line)[2])
        integral.append(sum/npts_integral)
vol=4.*math.pi*radius**3/4.
vol=cube_side**3
for iatom in range(num_atoms):
        print "For atom",iatom,"magnetic moment",integral[iatom]*vol






