shared = "../libDPG"
geometry = magnet.geo       
mesh     = magnet.vol.gz
constant geometryorder = 3

coefficient F (0,10,0), (0, 0, 0)

define fespace v -type=hcurlho_periodic -order=3
                 -yends=[0,1] -xends=[0,1] -dirichlet=[1]

define gridfunction u -fespace=v 

define bilinearform a -fespace=v -symmetric -spd
curlcurledge (1.0)
massedge  (0.001)


define linearform f -fespace=v
curledge F

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -solver=direct

define bilinearform acurl -fespace=v -symmetric -nonassemble
curlcurledge (1.0)

numproc drawflux np2 -bilinearform=acurl -solution=u  -label=flux

numproc visualization npv1 -vectorfunction=flux -clipsolution=vector -subdivision=3  -clipvec=[0,0,-1] 
