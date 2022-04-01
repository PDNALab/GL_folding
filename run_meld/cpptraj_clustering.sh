#! /bin/bash
#SBATCH --job-name=c      # Job name                            
#SBATCH --ntasks=1                # Number of MPI ranks               
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank   
#SBATCH --nodes=1                    # Number of nodes                
#SBATCH --ntasks-per-node=1                                           
#SBATCH --ntasks-per-socket=1                                         
#SBATCH --mem-per-cpu=50gb          # Memory per processor            
#SBATCH --output=hie.log
#SBATCH --time=12:00:00              # Time limit hrs:min:sec         
pwd; hostname; date                                                   
                                                                      
source ~/.load_Amber                                                  

for r in 2
do


    #Select residues on which to do the clustering
cat<<EOF>get_res.py
#! /usr/bin/env python

f = open('ss_nat.dat','r')

l = f.readlines()[0].rstrip()

ok = False
res = []
for i,s in enumerate(l):
    if ok and s not in 'HE':
        ok = False
        end = i
        res.append("{}-{}".format(ini,end))

    if (s in 'HE'):
        if not ok:
            ok = True
            ini = i + 1
        if (ok):
            pass


print ":{}".format(",".join(res))
EOF

res=`python ./get_res.py`
echo $res

#generate topology

rm -r Cpptraj_all_linkage_sieve_eps_$r
mkdir Cpptraj_all_linkage_sieve_eps_$r
cd Cpptraj_all_linkage_sieve_eps_$r

cat<<EFO>cpptraj.in
trajin ../trajectory.00.dcd 5000 300000 10
trajin ../trajectory.01.dcd 5000 300000 10
trajin ../trajectory.02.dcd 5000 300000 10
trajin ../trajectory.03.dcd 5000 300000 10
trajin ../trajectory.04.dcd 5000 300000 10
rms PredSSE first ${res}@CA,CB out trajrmsd.dat
cluster hieragglo epsilon $r linkage rms ${res}@CA,CB sieve 100 summary summary singlerepout representative repout unique repfmt pdb clusterout clusttraj avgout avg avgfmt pdb
go
EFO

for j in 0 1 2 3 4 5 6 7 8 9 
do
   cat<<EOF>>cpptraj.in
   clear trajin
   trajin  clusttraj.c$j
   trajout clusttraj.c$j.pdb model
   average average.c$j.pdb pdb
   reference unique.c$j.pdb [uniq$j]
   rms centr$j ref [uniq$j] @CA 
   atomicfluct out back.$j.apf @C,CA,N,O byres
   go
EOF
done


cpptraj -p ../../../system.top -i cpptraj.in
cd ..
done


