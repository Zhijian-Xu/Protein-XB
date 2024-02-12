from pymol import cmd
from pymol import stored
import pandas as pd
import com
import math

def calcDistance(lx,ly,lz,px,py,pz):
  return math.sqrt((lx-px)**2+(ly-py)**2+(lz-pz)**2);

def calcVectorLength(*vector):
  len = math.sqrt((vector[0])**2+(vector[1])**2+(vector[2])**2);
  return len

def multiply2vectors(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  return (vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2])

def calcVectorsAngle(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  len1=calcVectorLength(*vector1)
  len2=calcVectorLength(*vector2)
  if(len1 < 0.005):
     len1 = 0.005
  if(len2 < 0.005):
     len2 = 0.005
  dotMetrix = multiply2vectors(*vector_two)
  cos_angle = dotMetrix/(len1*len2)
  if(cos_angle > 1.0):
    cos_angle = 1.0
  angle = math.acos(cos_angle)
  return angle

def calcNormalVector(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  nv=[]
  nvx=vector1[1]*vector2[2]-vector1[2]*vector2[1]
  nvy=vector1[2]*vector2[0]-vector1[0]*vector2[2]
  nvz=vector1[0]*vector2[1]-vector1[1]*vector2[0]
  nv=[nvx,nvy,nvz]
  return nv

def normalize(*vector):
  vlen=calcVectorLength(*vector)
  if(vlen<0.005):
    vlen=0.005
  nvector=(vector[0]/vlen,vector[1]/vlen,vector[2]/vlen)
  return nvector

#get the coordinates of the selection
def get_coord_zj(selection_tmp):
  return cmd.get_coords(selection_tmp)[0]

def X_pi_distance_angle(sele_X,sele_CG,sele_D):
    sele_ring='(byres {}) and not name c+n+o+ca+cb+oh and not e. H'.format(sele_CG)
    com.com(sele_ring, object="sele2")
    sele2_xyz = get_coord_zj("sele2")
    CG_xyz = get_coord_zj(sele_ring + " and name CG")
    CD2_xyz = get_coord_zj(sele_ring + " and name CD2")
    vec1 = [sele2_xyz[0]-CG_xyz[0], sele2_xyz[1]-CG_xyz[1], sele2_xyz[2]-CG_xyz[2]]
    vec2 = [sele2_xyz[0]-CD2_xyz[0], sele2_xyz[1]-CD2_xyz[1], sele2_xyz[2]-CD2_xyz[2]]
    vec1_2 = vec1 + vec2
    nv=calcNormalVector(*vec1_2)
    
    sele1_xyz=get_coord_zj(sele_X)
    sele1_sele2 = [sele1_xyz[0]-sele2_xyz[0], sele1_xyz[1]-sele2_xyz[1], sele1_xyz[2]-sele2_xyz[2]] 
    alpha_vector = nv + sele1_sele2
    alpha = calcVectorsAngle(*alpha_vector)
    alpha = alpha*RAD_TO_DEG
    if(alpha > 90.0):
        alpha = 180.0 - alpha
    dst = cmd.distance("tmp", sele_X, "sele2")
    theta_ang = cmd.angle("tmp_angle", "sele2", sele_X, sele_D)
    
    return dst,alpha,theta_ang

def find_intra_halogenBond(PDB,halogen_chain,halogen_resi,halogen_name,cf):
    cmd.load(pdb_path+PDB+'.cif.gz')
    cmd.remove('solvent')
    halogen_bond_info=[]
    halogen_pairs=cmd.find_pairs('polymer.protein and c. {} and resi {} and name {}'.format(halogen_chain,halogen_resi,halogen_name),'polymer.protein and c. {} and (e. N+O+S or (r. phe+tyr+his+trp and name CG))'.format(halogen_chain),cutoff=cf)
    if halogen_pairs:
        for halogen_atom,pair_atom in halogen_pairs:
            cmd.select('sele_X','index {}'.format(halogen_atom[1]))
            cmd.select('sele_Y','index {}'.format(pair_atom[1]))
            cmd.select('sele_D','nbr. sele_x')

           
            atom_X=cmd.get_model('sele_X')
            for atom in atom_X.atom:
                chain_X=atom.chain
                resn_X=atom.resn
                resi_X=atom.resi
                name_X=atom.name
                symbol_X=atom.symbol
            index_X=str(halogen_atom[1])

            # print(symbol_X)
            atom_Y=cmd.get_model('sele_Y')
            for atom in atom_Y.atom:
                chain_Y=atom.chain
                resn_Y=atom.resn
                resi_Y=atom.resi
                name_Y=atom.name
                symbol_Y=atom.symbol
            index_Y=str(pair_atom[1])

            if name_Y=='CG':
                if symbol_X.upper()=='F':
                    max_dis=F_pi_dis
                elif symbol_X.upper()=='CL':
                    max_dis=Cl_pi_dis
                elif symbol_X.upper()=='BR':
                    max_dis=Br_pi_dis
                elif symbol_X.upper()=='I':
                    max_dis=I_pi_dis
                max_dis=max_dis+1 #extend
                halogen_dis,halogen_alphaang,halogen_thetaang=X_pi_distance_angle('sele_X','sele_Y','sele_D')

                if (min_dis<=halogen_dis<=max_dis) and (halogen_alphaang<X_pi_alpha) and (halogen_thetaang>X_pi_theta):
                    halogen_bond_info.append(((PDB,chain_X,resn_X,resi_X,name_X,symbol_X,index_X),(chain_Y,resn_Y,resi_Y,name_Y,'Pi',index_Y,str(max_dis-1),str(halogen_dis),str(halogen_alphaang)+'/'+str(halogen_thetaang))))
            

            else:

                halogen_dis=cmd.distance(selection1='sele_X',selection2='sele_Y')
                halogen_ang=cmd.angle(selection1='sele_D',selection2='sele_X',selection3='sele_Y')

                if symbol_X.upper()=='F':
                    r_X=van_rad_F
                elif symbol_X.upper()=='CL':
                    r_X=van_rad_Cl
                elif symbol_X.upper()=='BR':
                    r_X=van_rad_Br
                elif symbol_X.upper()=='I':
                    r_X=van_rad_I

                if symbol_Y.upper()=='C':
                    r_Y=van_rad_C
                elif symbol_Y.upper()=='N':
                    r_Y=van_rad_N
                elif symbol_Y.upper()=='O':
                    r_Y=van_rad_O
                elif symbol_Y.upper()=='S':
                    r_Y=van_rad_S
                
                max_dis=r_X+r_Y+1

                if (min_dis<=halogen_dis<=max_dis) and (min_angle<=halogen_ang<=max_angle):
                    halogen_bond_info.append(((PDB,chain_X,resn_X,resi_X,name_X,symbol_X,index_X),(chain_Y,resn_Y,resi_Y,name_Y,symbol_Y,index_Y,str(max_dis-1),str(halogen_dis),str(halogen_ang))))
    
    cmd.delete('all')
    return halogen_bond_info



materials_path='/home/dddc/jintian/protein_halogen/dataset/materials/'
pdb_path='/home/dddc/jintian/protein_halogen/dataset/data/'

van_rad_F=1.47
van_rad_Cl=1.75
van_rad_Br=1.85
van_rad_I=1.98

van_rad_C=1.70
van_rad_N=1.55
van_rad_O=1.52
van_rad_S=1.80

F_pi_dis=4.0
Cl_pi_dis=4.2
Br_pi_dis=4.3
I_pi_dis=4.5

X_pi_alpha=60+20
X_pi_theta=120-20


min_dis=0
cutoff=7.0
min_angle=120
max_angle=180

DEG_TO_RAD=3.1415926/180.0
RAD_TO_DEG=180.0/3.1415926


with open(materials_path+'protein_peptide_halogen_20220909_2433X_424PDB.txt','r')as f:
    all_halogen=f.readlines()

df_dic={'PDBid':[],'halogen_chain':[],'halogen_resn':[],'halogen_resi':[],'halogen_name':[],'halogen_symbol':[]}
for line in all_halogen:
    X_info=line.strip().split()
    df_dic['PDBid'].append(X_info[0])
    df_dic['halogen_chain'].append(X_info[1].split(':')[1])
    df_dic['halogen_resn'].append(X_info[3].split(':')[1])
    df_dic['halogen_resi'].append(X_info[2].split(':')[1])
    df_dic['halogen_name'].append(X_info[4].split(':')[1])
    df_dic['halogen_symbol'].append(X_info[5].split(':')[1])


halogen_info=pd.DataFrame(df_dic)
halogen_info=halogen_info.drop(index=halogen_info[halogen_info['PDBid']=='2ag3'].index.tolist())

input_param_ls=list(zip(halogen_info['PDBid'],halogen_info['halogen_chain'],halogen_info['halogen_resi'],halogen_info['halogen_name'],[cutoff]*halogen_info.shape[0]))


from multiprocessing import Pool
def multi_find_intra_halogenBond(param_ls):
    return find_intra_halogenBond(param_ls[0],param_ls[1],param_ls[2],param_ls[3],param_ls[4])

all_intra_halogenbond=[]
pool=Pool(128)
intra_halogenbonds=pool.map(multi_find_intra_halogenBond,input_param_ls)
pool.close()
pool.join()
for info in intra_halogenbonds:
    all_intra_halogenbond.extend(info)


with open(materials_path+'halogen_bond_chainsinfo2+1_231208_pi80-100.csv','w')as f:
    f.write(','.join(['PDBid','halogen_chain','halogen_resn','halogen_resi','halogen_name','halogen_symbol','halogen_index','pair_chain','pair_resn','pair_resi','pair_name','pair_symbol','pair_index','vander_R','distance','angle'])+'\n')
    for info_X,info_Y in all_intra_halogenbond:
        f.write(','.join(info_X+info_Y)+'\n')

