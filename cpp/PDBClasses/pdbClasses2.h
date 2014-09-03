#ifndef _PDBCLASSES_
#define _PDBCLASSES_

using namespace std;
#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#define PI 3.141592653589793

#define DBL (10e10)
/* _CA_ alpha carbon atom
   _CB_ beta carbon atom
   _HV_ heavy atom
   _HN_ backbone hydrogen bond
   _METHYLC_ carbon in a methyl group
*/
typedef enum { _CA_,_CB_,_HV_,_HN_,_METHYLC_} gen_atom_type;
/* _CA_CA_ contacts between CA atoms    _CB_ contacts between CB atoms
  _HV_ contacts between any heavy atom    
  _HVSC_HVSC_ contacts between heavy atoms in the sidechain
  _HVBK_HVBK_ contacts between heavy atoms in the backbone
  _HVBK_HVSC_ contact between one heavy atom in the backbone, the other
            in the sidechain. Order does not matter.
  _F_HVBK_S_HVSC_ contact between one heavy atom in the backbone of the
first amino acid  with one heavy atom in the sidechain of the second 
amino acid
  _F_HVSC_S_HVBK_  contact between one heavy atom in the sidechain of the
first amino acid  with one heavy atom in the backbone of the second
amino acid
  _HN_HN_ contact between hydrogens in the backbone. Return 0 if the
H are not explicitly represented.
  _HN_METHYLC_ contact between a hydrogen in the backbone and a carbon belonging to a methyl group
  _HN_METHYLC_ the reverse of _HN_METHYLC_
  _METHYLC_METHYLC_ contact between carbons belonging to two methyl groups
*/
typedef enum {_CA_CA_=0,_CB_CB_,_HV_HV_,_HVSC_HVSC_,_HVBK_HVBK_,_HVBK_HVSC_,_F_HVBK_S_HVSC_,_F_HVSC_S_HVBK_,_HN_HN_,_HN_METHYLC_,_METHYLC_HN_,_METHYLC_METHYLC_} cont_sel; /*method to select a contact*/
typedef enum {_ALPHA_=0, _BETA_, _ALPHA_PLUS_BETA_, _ALPHA_BETA_, _MULTI_DOM,
	      _MEMBR_CELL_SF_, _SMALL_} fold; 
typedef enum {_NONE_=0,_Z_,_gz_,_bz2_,_zip_} zip_format;
static const string all_aa="ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL";
static const int Nside_chain_atoms = 216 ;
static string side_chain_atoms[ Nside_chain_atoms ]={
"ATOM     57  CB  ALA C  12      -0.889   0.050   1.245  1.00  0.34           C  ",
"ATOM     59  HA  ALA C  12      -0.620   0.029  -0.883  1.00  0.32           H  ",
"ATOM     60 1HB  ALA C  12      -0.309  -0.226   2.113  1.00  1.11           H  ",
"ATOM     61 2HB  ALA C  12      -1.272   1.053   1.370  1.00  1.05           H  ",
"ATOM     62 3HB  ALA C  12      -1.712  -0.639   1.126  1.00  1.06           H  ",
"ATOM    166  CB  ARG C  19      -0.878   0.042   1.254  1.00  0.70           C  ",
"ATOM    167  CG  ARG C  19      -2.002  -0.990   1.128  1.00  0.95           C  ",
"ATOM    168  CD  ARG C  19      -3.353  -0.305   1.342  1.00  1.11           C  ",
"ATOM    169  NE  ARG C  19      -4.347  -0.853   0.377  1.00  1.32           N  ",
"ATOM    170  CZ  ARG C  19      -5.476  -0.229   0.173  1.00  1.79           C  ",
"ATOM    171  NH1 ARG C  19      -5.736   0.881   0.812  1.00  2.38           N  ",
"ATOM    172  NH2 ARG C  19      -6.345  -0.714  -0.671  1.00  2.22           N  ",
"ATOM    174  HA  ARG C  19      -0.623   0.029  -0.880  1.00  0.58           H  ",
"ATOM    175 1HB  ARG C  19      -0.276  -0.185   2.122  1.00  0.75           H  ",
"ATOM    176 2HB  ARG C  19      -1.307   1.027   1.360  1.00  0.86           H  ",
"ATOM    177 1HG  ARG C  19      -1.976  -1.435   0.144  1.00  1.49           H  ",
"ATOM    178 2HG  ARG C  19      -1.870  -1.758   1.874  1.00  1.53           H  ",
"ATOM    179 1HD  ARG C  19      -3.693  -0.488   2.351  1.00  1.68           H  ",
"ATOM    180 2HD  ARG C  19      -3.246   0.758   1.184  1.00  1.52           H  ",
"ATOM    181  HE  ARG C  19      -4.154  -1.685  -0.104  1.00  1.68           H  ",
"ATOM    182 1HH1 ARG C  19      -5.072   1.255   1.458  1.00  2.38           H  ",
"ATOM    183 2HH1 ARG C  19      -6.601   1.357   0.653  1.00  3.05           H  ",
"ATOM    184 1HH2 ARG C  19      -6.146  -1.563  -1.161  1.00  2.32           H  ",
"ATOM    185 2HH2 ARG C  19      -7.209  -0.237  -0.828  1.00  2.72           H  ",
"ATOM    794  CB  ASN C  59      -0.885   0.029   1.243  1.00  0.39           C  ",
"ATOM    795  CG  ASN C  59      -0.013   0.213   2.480  1.00  0.38           C  ",
"ATOM    796  OD1 ASN C  59       1.149  -0.134   2.469  1.00  0.54           O  ",
"ATOM    797  ND2 ASN C  59      -0.529   0.745   3.554  1.00  0.42           N  ",
"ATOM    799  HA  ASN C  59      -0.621   0.029  -0.883  1.00  0.34           H  ",
"ATOM    800 1HB  ASN C  59      -1.580   0.847   1.167  1.00  0.45           H  ",
"ATOM    801 2HB  ASN C  59      -1.428  -0.900   1.321  1.00  0.45           H  ",
"ATOM    802 1HD2 ASN C  59       0.024   0.872   4.352  1.00  0.43           H  ",
"ATOM    803 2HD2 ASN C  59      -1.469   1.020   3.561  1.00  0.56           H  ",
"ATOM    107  CB  ASP C  15      -0.877   0.053   1.257  1.00  0.53           C  ",
"ATOM    108  CG  ASP C  15      -1.612  -1.278   1.443  1.00  0.90           C  ",
"ATOM    109  OD1 ASP C  15      -0.951  -2.267   1.713  1.00  1.79           O  ",
"ATOM    110  OD2 ASP C  15      -2.825  -1.283   1.317  1.00  1.17           O  ",
"ATOM    112  HA  ASP C  15      -0.618   0.027  -0.884  1.00  0.48           H  ",
"ATOM    113 1HB  ASP C  15      -0.256   0.244   2.120  1.00  0.82           H  ",
"ATOM    114 2HB  ASP C  15      -1.601   0.848   1.155  1.00  0.88           H  ",
"ATOM    255  CB  CYS    16      -0.979  -0.052   1.171  1.00  3.44           C  ",
"ATOM    256  SG  CYS    16      -1.946   1.449   1.431  1.00  2.74           S  ",
"ATOM    258  HA  CYS    16      -0.562   0.047  -0.797  1.00  1.34           H  ",
"ATOM    259 1HB  CYS    16      -1.621  -0.776   1.027  1.00  1.04           H  ",
"ATOM    260 2HB  CYS    16      -0.503  -0.230   2.008  1.00  1.00           H  ",
"ATOM    691  CB  GLN C  52      -0.881   0.040   1.252  1.00  1.16           C  ",
"ATOM    692  CG  GLN C  52      -1.990  -1.008   1.133  1.00  1.63           C  ",
"ATOM    693  CD  GLN C  52      -3.101  -0.692   2.136  1.00  2.38           C  ",
"ATOM    694  OE1 GLN C  52      -2.839  -0.480   3.304  1.00  2.78           O  ",
"ATOM    695  NE2 GLN C  52      -4.339  -0.651   1.728  1.00  3.19           N  ",
"ATOM    697  HA  GLN C  52      -0.624   0.029  -0.881  1.00  1.01           H  ",
"ATOM    698 1HB  GLN C  52      -0.277  -0.170   2.123  1.00  1.26           H  ",
"ATOM    699 2HB  GLN C  52      -1.323   1.021   1.348  1.00  1.17           H  ",
"ATOM    700 1HG  GLN C  52      -2.394  -0.992   0.130  1.00  1.89           H  ",
"ATOM    701 2HG  GLN C  52      -1.585  -1.986   1.343  1.00  1.83           H  ",
"ATOM    702 1HE2 GLN C  52      -5.059  -0.450   2.364  1.00  3.80           H  ",
"ATOM    703 2HE2 GLN C  52      -4.551  -0.821   0.787  1.00  3.44           H  ",
"ATOM    498  CB  GLU C  39      -0.882   0.041   1.251  1.00  0.60           C  ",
"ATOM    499  CG  GLU C  39      -1.861   1.214   1.150  1.00  0.99           C  ",
"ATOM    500  CD  GLU C  39      -2.952   0.884   0.129  1.00  1.86           C  ",
"ATOM    501  OE1 GLU C  39      -3.707  -0.042   0.377  1.00  2.53           O  ",
"ATOM    502  OE2 GLU C  39      -3.014   1.564  -0.882  1.00  2.58           O  ",
"ATOM    504  HA  GLU C  39      -0.624   0.029  -0.881  1.00  0.62           H  ",
"ATOM    505 1HB  GLU C  39      -1.435  -0.883   1.331  1.00  1.01           H  ",
"ATOM    506 2HB  GLU C  39      -0.262   0.167   2.125  1.00  0.93           H  ",
"ATOM    507 1HG  GLU C  39      -2.312   1.390   2.116  1.00  1.56           H  ",
"ATOM    508 2HG  GLU C  39      -1.330   2.099   0.834  1.00  1.55           H  ",
"ATOM    323 1HA  GLY    20      -0.571   0.038  -0.800  1.00  1.87           H  ",
"ATOM    324 2HA  GLY    20      -0.620   0.027   0.781  1.00  1.04           H  ",
"ATOM    609  CB  HIS C  46      -0.885   0.033   1.250  1.00  0.43           C  ",
"ATOM    610  CG  HIS C  46      -1.431   1.420   1.451  1.00  0.55           C  ",
"ATOM    611  ND1 HIS C  46      -2.473   1.922   0.688  1.00  0.71           N  ",
"ATOM    612  CD2 HIS C  46      -1.088   2.423   2.323  1.00  0.77           C  ",
"ATOM    613  CE1 HIS C  46      -2.718   3.175   1.111  1.00  0.79           C  ",
"ATOM    614  NE2 HIS C  46      -1.902   3.530   2.107  1.00  0.84           N  ",
"ATOM    616  HA  HIS C  46      -0.623   0.030  -0.882  1.00  0.39           H  ",
"ATOM    617 1HB  HIS C  46      -1.704  -0.662   1.129  1.00  0.48           H  ",
"ATOM    618 2HB  HIS C  46      -0.299  -0.250   2.112  1.00  0.46           H  ",
"ATOM    619  HD1 HIS C  46      -2.946   1.450  -0.030  1.00  0.87           H  ",
"ATOM    620  HD2 HIS C  46      -0.305   2.362   3.064  1.00  0.98           H  ",
"ATOM    621  HE1 HIS C  46      -3.481   3.816   0.695  1.00  0.95           H  ",
"ATOM    421  CB  ILE C  34      -0.886   0.040   1.247  1.00  0.33           C  ",
"ATOM    422  CG1 ILE C  34      -2.009  -0.990   1.105  1.00  0.39           C  ",
"ATOM    423  CG2 ILE C  34      -1.492   1.438   1.402  1.00  0.36           C  ",
"ATOM    424  CD1 ILE C  34      -2.967  -0.554  -0.006  1.00  0.41           C  ",
"ATOM    426  HA  ILE C  34      -0.621   0.028  -0.883  1.00  0.31           H  ",
"ATOM    427  HB  ILE C  34      -0.289  -0.191   2.119  1.00  0.37           H  ",
"ATOM    428 1HG1 ILE C  34      -1.584  -1.953   0.859  1.00  0.42           H  ",
"ATOM    429 2HG1 ILE C  34      -2.551  -1.063   2.037  1.00  0.43           H  ",
"ATOM    430 1HG2 ILE C  34      -0.709   2.179   1.346  1.00  1.07           H  ",
"ATOM    431 2HG2 ILE C  34      -2.208   1.610   0.613  1.00  0.99           H  ",
"ATOM    432 3HG2 ILE C  34      -1.988   1.511   2.359  1.00  1.03           H  ",
"ATOM    433 1HD1 ILE C  34      -3.317   0.449   0.194  1.00  1.06           H  ",
"ATOM    434 2HD1 ILE C  34      -2.451  -0.573  -0.954  1.00  1.16           H  ",
"ATOM    435 3HD1 ILE C  34      -3.810  -1.229  -0.041  1.00  1.08           H  ",
"ATOM     67  CB  LEU C  13      -0.869   0.042   1.259  1.00  0.37           C  ",
"ATOM     68  CG  LEU C  13      -2.069  -0.889   1.089  1.00  0.43           C  ",
"ATOM     69  CD1 LEU C  13      -2.775  -1.063   2.435  1.00  0.52           C  ",
"ATOM     70  CD2 LEU C  13      -3.043  -0.281   0.077  1.00  0.54           C  ",
"ATOM     72  HA  LEU C  13      -0.631   0.036  -0.874  1.00  0.39           H  ",
"ATOM     73 1HB  LEU C  13      -0.283  -0.276   2.110  1.00  0.41           H  ",
"ATOM     74 2HB  LEU C  13      -1.218   1.050   1.422  1.00  0.44           H  ",
"ATOM     75  HG  LEU C  13      -1.731  -1.852   0.733  1.00  0.48           H  ",
"ATOM     76 1HD1 LEU C  13      -2.652  -0.165   3.025  1.00  1.12           H  ",
"ATOM     77 2HD1 LEU C  13      -3.826  -1.244   2.270  1.00  1.15           H  ",
"ATOM     78 3HD1 LEU C  13      -2.342  -1.901   2.962  1.00  1.20           H  ",
"ATOM     79 1HD2 LEU C  13      -2.943   0.794   0.084  1.00  1.27           H  ",
"ATOM     80 2HD2 LEU C  13      -2.820  -0.658  -0.909  1.00  1.17           H  ",
"ATOM     81 3HD2 LEU C  13      -4.055  -0.550   0.345  1.00  1.06           H  ",
"ATOM    295  CB  LYS C  27      -0.882   0.043   1.249  1.00  0.68           C  ",
"ATOM    296  CG  LYS C  27      -1.941  -1.057   1.166  1.00  1.30           C  ",
"ATOM    297  CD  LYS C  27      -2.584  -1.248   2.540  1.00  1.84           C  ",
"ATOM    298  CE  LYS C  27      -4.022  -0.727   2.507  1.00  2.70           C  ",
"ATOM    299  NZ  LYS C  27      -4.717  -1.106   3.770  1.00  3.26           N  ",
"ATOM    301  HA  LYS C  27      -0.624   0.027  -0.880  1.00  0.65           H  ",
"ATOM    302 1HB  LYS C  27      -0.271  -0.109   2.127  1.00  1.13           H  ",
"ATOM    303 2HB  LYS C  27      -1.369   1.005   1.312  1.00  0.94           H  ",
"ATOM    304 1HG  LYS C  27      -2.698  -0.775   0.449  1.00  1.80           H  ",
"ATOM    305 2HG  LYS C  27      -1.477  -1.981   0.857  1.00  1.89           H  ",
"ATOM    306 1HD  LYS C  27      -2.587  -2.299   2.794  1.00  2.25           H  ",
"ATOM    307 2HD  LYS C  27      -2.022  -0.699   3.281  1.00  2.16           H  ",
"ATOM    308 1HE  LYS C  27      -4.013   0.348   2.409  1.00  3.16           H  ",
"ATOM    309 2HE  LYS C  27      -4.542  -1.161   1.666  1.00  3.07           H  ",
"ATOM    310 1HZ  LYS C  27      -4.112  -0.873   4.582  1.00  3.62           H  ",
"ATOM    311 2HZ  LYS C  27      -5.613  -0.583   3.841  1.00  3.55           H  ",
"ATOM    312 3HZ  LYS C  27      -4.911  -2.128   3.766  1.00  3.54           H  ",
"ATOM    491  CB  MET A 181      -0.905   0.057   1.230  1.00  5.55      1CKA 570",
"ATOM    492  CG  MET A 181      -1.874  -1.109   1.299  1.00  7.39      1CKA 571",
"ATOM    493  SD  MET A 181      -2.865  -1.296  -0.220  1.00 14.41      1CKA 572",
"ATOM    494  CE  MET A 181      -3.803   0.232  -0.118  1.00 19.39      1CKA 573",
"ATOM     21  CB  PHE C  10      -0.880   0.038   1.253  1.00  0.34           C  ",
"ATOM     22  CG  PHE C  10      -2.152  -0.739   1.005  1.00  0.35           C  ",
"ATOM     23  CD1 PHE C  10      -2.094  -2.026   0.457  1.00  1.20           C  ",
"ATOM     24  CD2 PHE C  10      -3.392  -0.170   1.322  1.00  1.31           C  ",
"ATOM     25  CE1 PHE C  10      -3.274  -2.744   0.227  1.00  1.21           C  ",
"ATOM     26  CE2 PHE C  10      -4.572  -0.888   1.090  1.00  1.34           C  ",
"ATOM     27  CZ  PHE C  10      -4.513  -2.174   0.542  1.00  0.45           C  ",
"ATOM     29  HA  PHE C  10      -0.625   0.029  -0.880  1.00  0.32           H  ",
"ATOM     30 1HB  PHE C  10      -0.345  -0.404   2.081  1.00  0.41           H  ",
"ATOM     31 2HB  PHE C  10      -1.125   1.063   1.489  1.00  0.33           H  ",
"ATOM     32  HD1 PHE C  10      -1.139  -2.467   0.215  1.00  2.11           H  ",
"ATOM     33  HD2 PHE C  10      -3.439   0.821   1.745  1.00  2.21           H  ",
"ATOM     34  HE1 PHE C  10      -3.229  -3.737  -0.196  1.00  2.10           H  ",
"ATOM     35  HE2 PHE C  10      -5.528  -0.447   1.334  1.00  2.25           H  ",
"ATOM     36  HZ  PHE C  10      -5.423  -2.727   0.364  1.00  0.50           H  ",
"ATOM    769  CB  PRO C  57      -0.814  -0.092   1.293  1.00  0.19           C  ",
"ATOM    770  CG  PRO C  57      -0.057  -1.073   2.218  1.00  0.21           C  ",
"ATOM    771  CD  PRO C  57       0.847  -1.927   1.308  1.00  0.25           C  ",
"ATOM    772  HA  PRO C  57      -0.660   0.030  -0.852  1.00  0.22           H  ",
"ATOM    773 1HB  PRO C  57      -0.884   0.880   1.754  1.00  0.20           H  ",
"ATOM    774 2HB  PRO C  57      -1.800  -0.476   1.084  1.00  0.19           H  ",
"ATOM    775 1HG  PRO C  57       0.540  -0.520   2.931  1.00  0.24           H  ",
"ATOM    776 2HG  PRO C  57      -0.760  -1.709   2.733  1.00  0.22           H  ",
"ATOM    777 1HD  PRO C  57       1.859  -1.943   1.690  1.00  0.29           H  ",
"ATOM    778 2HD  PRO C  57       0.458  -2.932   1.220  1.00  0.26           H  ",
"ATOM    626  CB  SER C  47      -0.882   0.041   1.249  1.00  0.94           C  ",
"ATOM    627  OG  SER C  47      -1.817  -1.028   1.199  1.00  1.13           O  ",
"ATOM    629  HA  SER C  47      -0.624   0.030  -0.882  1.00  0.82           H  ",
"ATOM    630 1HB  SER C  47      -0.270  -0.065   2.129  1.00  0.90           H  ",
"ATOM    631 2HB  SER C  47      -1.404   0.989   1.289  1.00  1.20           H  ",
"ATOM    632  HG  SER C  47      -2.495  -0.860   1.858  1.00  1.44           H  ",
"ATOM    219  CB  THR C  22      -0.880   0.042   1.251  1.00  0.51           C  ",
"ATOM    220  OG1 THR C  22      -0.660  -1.134   2.019  1.00  0.64           O  ",
"ATOM    221  CG2 THR C  22      -2.351   0.118   0.841  1.00  0.70           C  ",
"ATOM    223  HA  THR C  22      -0.626   0.030  -0.879  1.00  0.49           H  ",
"ATOM    224  HB  THR C  22      -0.629   0.910   1.841  1.00  0.53           H  ",
"ATOM    225  HG1 THR C  22      -1.028  -0.992   2.893  1.00  1.07           H  ",
"ATOM    226 1HG2 THR C  22      -2.510   0.995   0.232  1.00  1.33           H  ",
"ATOM    227 2HG2 THR C  22      -2.615  -0.765   0.277  1.00  1.26           H  ",
"ATOM    228 3HG2 THR C  22      -2.968   0.177   1.726  1.00  1.23           H  ",
"ATOM    532  CB  TRP C  42      -0.892   0.034   1.244  1.00  0.30           C  ",
"ATOM    533  CG  TRP C  42      -1.995  -0.965   1.089  1.00  0.27           C  ",
"ATOM    534  CD1 TRP C  42      -1.859  -2.302   1.248  1.00  0.28           C  ",
"ATOM    535  CD2 TRP C  42      -3.391  -0.733   0.743  1.00  0.28           C  ",
"ATOM    536  NE1 TRP C  42      -3.084  -2.907   1.023  1.00  0.28           N  ",
"ATOM    537  CE2 TRP C  42      -4.060  -1.981   0.708  1.00  0.27           C  ",
"ATOM    538  CE3 TRP C  42      -4.137   0.426   0.461  1.00  0.33           C  ",
"ATOM    539  CZ2 TRP C  42      -5.418  -2.073   0.402  1.00  0.30           C  ",
"ATOM    540  CZ3 TRP C  42      -5.505   0.336   0.154  1.00  0.36           C  ",
"ATOM    541  CH2 TRP C  42      -6.143  -0.912   0.124  1.00  0.34           C  ",
"ATOM    543  HA  TRP C  42      -0.619   0.029  -0.885  1.00  0.24           H  ",
"ATOM    544 1HB  TRP C  42      -0.302  -0.209   2.116  1.00  0.34           H  ",
"ATOM    545 2HB  TRP C  42      -1.314   1.022   1.357  1.00  0.35           H  ",
"ATOM    546  HD1 TRP C  42      -0.944  -2.815   1.509  1.00  0.31           H  ",
"ATOM    547  HE1 TRP C  42      -3.258  -3.871   1.077  1.00  0.30           H  ",
"ATOM    548  HE3 TRP C  42      -3.654   1.393   0.481  1.00  0.36           H  ",
"ATOM    549  HZ2 TRP C  42      -5.906  -3.035   0.378  1.00  0.31           H  ",
"ATOM    550  HZ3 TRP C  42      -6.068   1.233  -0.062  1.00  0.43           H  ",
"ATOM    551  HH2 TRP C  42      -7.194  -0.976  -0.111  1.00  0.38           H  ",
"ATOM    119  CB  TYR C  16      -0.908   0.037   1.227  1.00  0.47           C  ",
"ATOM    120  CG  TYR C  16      -1.455   1.421   1.440  1.00  0.45           C  ",
"ATOM    121  CD1 TYR C  16      -2.248   2.000   0.455  1.00  1.32           C  ",
"ATOM    122  CD2 TYR C  16      -1.178   2.116   2.621  1.00  1.14           C  ",
"ATOM    123  CE1 TYR C  16      -2.772   3.284   0.640  1.00  1.34           C  ",
"ATOM    124  CE2 TYR C  16      -1.700   3.400   2.814  1.00  1.13           C  ",
"ATOM    125  CZ  TYR C  16      -2.498   3.986   1.821  1.00  0.45           C  ",
"ATOM    126  OH  TYR C  16      -3.011   5.253   2.010  1.00  0.48           O  ",
"ATOM    128  HA  TYR C  16      -0.610   0.018  -0.883  1.00  0.51           H  ",
"ATOM    129 1HB  TYR C  16      -1.736  -0.634   1.067  1.00  0.51           H  ",
"ATOM    130 2HB  TYR C  16      -0.353  -0.270   2.100  1.00  0.48           H  ",
"ATOM    131  HD1 TYR C  16      -2.457   1.449  -0.449  1.00  2.17           H  ",
"ATOM    132  HD2 TYR C  16      -0.564   1.661   3.383  1.00  1.99           H  ",
"ATOM    133  HE1 TYR C  16      -3.389   3.733  -0.126  1.00  2.21           H  ",
"ATOM    134  HE2 TYR C  16      -1.491   3.938   3.728  1.00  1.97           H  ",
"ATOM    135  HH  TYR C  16      -2.483   5.866   1.492  1.00  0.84           H  ",
"ATOM    239  CB  VAL    15      -1.055   0.022   1.109  1.00  2.71           C  ",
"ATOM    240  CG1 VAL    15      -1.650   1.393   1.276  1.00  4.19           C  ",
"ATOM    241  CG2 VAL    15      -2.155  -0.998   0.843  1.00  3.20           C  ",
"ATOM    243  HA  VAL    15      -0.524   0.022  -0.832  1.00  3.09           H  ",
"ATOM    244  HB  VAL    15      -0.573  -0.217   1.942  1.00  2.45           H  ",
"ATOM    245 1HG1 VAL    15      -1.317   2.026   0.605  1.00  7.84           H  ",
"ATOM    246 2HG1 VAL    15      -2.631   1.379   1.184  1.00  2.23           H  ",
"ATOM    247 3HG1 VAL    15      -1.465   1.785   2.160  1.00  2.26           H  ",
"ATOM    248 1HG2 VAL    15      -1.825  -1.898   0.677  1.00  7.04           H  ",
"ATOM    249 2HG2 VAL    15      -2.780  -1.044   1.591  1.00  1.00           H  ",
"ATOM    250 3HG2 VAL    15      -2.699  -0.731   0.076  1.00  1.52           H  "
};

/*define the GEOMETRY for each amino acid*/

const static double rpi = PI/180.0;

const static double N_CA_CB    = 110.4;
const static double CB_N_CA_C  = 122.0;
const static double N_CA_CBXY  = 180.0-acos(cos(N_CA_CB*rpi)/cos(CB_N_CA_C*rpi))/rpi;
const static double PRO_N_CA_CB  = 103.0;
const static double PRO_N_CA_CBXY= 180.0-acos(cos(PRO_N_CA_CB*rpi)/cos(CB_N_CA_C*rpi))/rpi;

const static double N_CA_C   = 111.2;
const static double CA_C_O   = 120.8;
const static double C_N_CA   = 121.7;
const static double CA_C_N   = 116.2;

const static double C_CA_CBXY = 360.0-N_CA_CBXY-N_CA_C;
const static double C_CA_CB   = 180.0-acos(cos(CB_N_CA_C*rpi)*cos(C_CA_CBXY*rpi))/rpi;
const static double PRO_C_CA_CBXY = 360.0-PRO_N_CA_CBXY-N_CA_C;
const static double PRO_C_CA_CB   = 180.0-acos(cos(CB_N_CA_C*rpi)*cos(PRO_C_CA_CBXY*rpi))/rpi;

const static double N_CA     = 1.458;
const static double CA_C     = 1.525;
const static double C_O      = 1.231;
const static double N_H      = 1.020;
const static double C_N      = 1.329;
const static double CA_CB    = 1.521;
/*Determine the intrinsic axis of the amino acids
  FOR INITIALIZING AMIN ACID STRUCTURE
  
  Two Streched amino acids(psi = -180, and phi = 180)
                     (ix)
       Ca--(iy)   N   |    C
     / a|  \    /  \ a|  /
    N   |    C       Ca--(iy)
       (ix)   
  angle : a is such that the two intrinsic x axis are parallel
  ------------------------------*/
const static double C_CA     = sqrt(C_N*C_N + N_CA*N_CA - 2.0*C_N*N_CA*cos(C_N_CA*rpi));
const static double N_C0_CA  = asin(N_CA*sin(C_N_CA*rpi)/C_CA)/rpi;
const static double C_CA_N0  = 180.0 - N_C0_CA - C_N_CA;
const static double CA_C_CA  = CA_C_N + N_C0_CA;
const static double CA_CA = sqrt(CA_C*CA_C + C_CA*C_CA - 2.0*CA_C*C_CA*cos(CA_C_CA*rpi));
const static double CA_CA_C0 = asin(CA_C*sin(CA_C_CA*rpi)/CA_CA)/rpi;
const static double C_CA0_CA = 180.0 - CA_C_CA - CA_CA_C0;
const static double CA_CA_N0 = C_CA_N0 - CA_CA_C0;

const static double init_N_CA_x   = (N_CA_C + C_CA0_CA - CA_CA_N0)/2.0;
const static double init_C_CA_x   = N_CA_C - init_N_CA_x;

const static double init_N_CA_y   = init_N_CA_x + 90 ;
const static double init_C_CA_y   = 90 - init_C_CA_x;
const static double init_y_C_CA   = 180 - init_C_CA_y;
const static double init_O_C_y    = 360 - init_y_C_CA - CA_C_O;
/*===========================================================*/
/*Functions that deal with generalized atom types*/
string gat2str(gen_atom_type );
gen_atom_type str2gat(string &);
/*===========================================================*/
/*Functions that deal with contact selections*/
string cs2str(cont_sel &);
cont_sel str2cs(string &);
bool symmetric_map(cont_sel &);
/*===========================================================*/
/*Functions to print maps*/
void print_symmetric_map(ostream &, double **, int &);
void print_asymmetric_map(ostream &, double **, int &);
/*===========================================================*/
/*Functions that deal with zipped formats*/
zip_format returnZipFormat(const char *);
string zf2str(zip_format);
bool unzip(char *,const char *);
bool zip(char *,const char *,zip_format &);
bool zip(const char *,zip_format &);
/*===========================================================*/
double det_3x3( double ** ) ;
/*===========================================================*/
void remove_bound_cond(double **,const int&,const double&);
/*===========================================================*/
/*Some funtions to convert strings to numbers and viceversa*/
void int2string5( const int &, string & ) ;
void int2string4( const int &, string & ) ;
void string2int( const string &,  int & ) ;
void double2string8_3( const double & , string & ) ;
void double2string6_2( const double & , string & ) ;
void string2double( const string &, double & ) ;
/*----some handy functions-------------------*/
bool is_ATOM( string & ) ;
bool is_ANISOU( string & ) ;
bool is_TER( string & ) ;
bool is_END( string & ) ;
bool is_END( ifstream & ) ;
bool is_HETATM( string & ) ;
bool is_SIGATM( string & ) ;
bool is_AMINOACID( string & ) ;
string three_letter_code_to_one_letter_code( const string &) ;
void dump_array( double **, int, int ) ;
void dump_array2( double **, int, int ) ;
void bubble_sort( double [ ], const int & ) ;
void bubble_sort( double [ ], const int &, double [ ] ) ;
double **alloc_array( int, int ) ;
void init_array(const double &,double **,int &);
void assign_array(double **,double **,int &);
void mult_arrayBA(double **,double **,int &);
void mult_arrayCBA(double **,double **,double **,int &);
int **alloc_int_array( int, int ) ;
void init_int_array(int &,int **,int &);
void assign_int_array(int **,int **,int &);
void mult_int_arrayBA(int **,int **,int &);
void mult_int_arrayCBA(int **,int **,int **,int &);
double get_rot_matrix_by_rms( double **, double **, long, double ** ) ;
double get_rmsII(double **, double **, double **, long ) ;
double get_rms(double **, double **, long ) ;
/*======================================================*/
class PDBvector{

  friend string  &operator<<(string  &, const PDBvector &);//with pdb format
  friend ostream &operator<<(ostream &, const PDBvector &);//with pdb format

 public:
  PDBvector( ) ;
  PDBvector(double const &, double const & , double const &) ;
  PDBvector( double * ) ;
  void assign( double const &, double const & , double const &) ;
  void assign( double * ) ;
  double &operator[]( const int & );  
  double getComponent( const int & );
  double *getComponents( ) ;
  void getComponentsII( double * ) ;
  PDBvector &operator=( const PDBvector & ) ;
  PDBvector &operator+=( const PDBvector & ) ;
  PDBvector operator+( const PDBvector & ) ;
  PDBvector &operator-=( const PDBvector & ) ;
  PDBvector operator-( const PDBvector & ) ;
  PDBvector &operator/=( const double & ) ;
  PDBvector operator/( const double & ) ;
  PDBvector &operator*=( const double & ) ;
  PDBvector operator*( const double & ) ;
  double operator*( const PDBvector & ) ;
  double d2( const PDBvector & ) ;
  double norm2( ) ;
  void normalize( ) ;
  double angle_with( PDBvector & ) ;
  PDBvector operator^( const PDBvector & ) ;
  void rotate( const PDBvector &, const double &);/*rotate around axis some degrees*/
  PDBvector rotateII( const PDBvector &, const double &);
  void rotate( const PDBvector &, const PDBvector &, const double &);/*rotate at
                                             new origin around axis some degrees*/
  PDBvector rotateII( const PDBvector &, const PDBvector &, const double &);
  void rotate( double **  );
  PDBvector rotateII( double ** );
  double radius(PDBvector&, PDBvector&);/*radius of the circle that
goes through the three vectors*/

 protected:
  double x ;
  double y ;
  double z ;  

} ;


bool isPDBatom( string &) ;
//------------ class PDBatom  --------------


class PDBatom{            //same information than PDB files

  friend class listPDBatom ;
  friend class DMDAtom ;

  friend string &operator<<(string &, const PDBatom &);
  friend ostream &operator<<( ostream &, const PDBatom & ) ;
  friend string &operator>>(string &,       PDBatom &);
  friend ifstream &operator>>(ifstream &, PDBatom &);

 public:
  PDBatom ( ) ;
  PDBatom( const PDBvector & ) ;
  PDBvector &getCoord( ) ;
  string exportAtom(int mode=0);
  bool is_recorded( ) ;
  bool is_backbone( ) ;
  bool matchresSeq( const PDBatom & ) ;
  void empty( ) ;
  string *printTER( ) ;
  bool isHeavy( ) ;
  bool doContact( const PDBatom & , const double & ) ;
  double d2( PDBatom & ) ;
  int getResSeq( ) ;
  string getiCode();
  void setiCode(const string &);
  string * printName( ) ;
  void getAtomName( string & ) ;
  string getAtomName();
  string * printResName( ) ;
  string getResName( ) ;
  string * printChainID( ) ;
  void changeResSeq( const int & ) ;
  void printChainID( string & ) ;
  void setChainID( string & ) ;
  void assignCoord( double * ) ;
  void assignCoord( PDBvector & ) ;
  bool is( const string & ) ;
  bool is_gen_atom_type( gen_atom_type );
  double *getCoordToDoubleArray( ) ;
  void getCoordToDoubleArrayII(double * ) ;
  void changeAtomSerialNumber( const int & ) ;
  int getSerialNumber( ) ;
  void assignName( string & ) ;
  void assignResName( string & ) ;
  double angle_with( PDBatom&, PDBatom& ) ;
  PDBvector PDBvector_to( PDBatom & ) ;
  void translate( const PDBvector & ) ;
  void rotate( const PDBvector &, const double & );/*rotate around axis some 
                                                     degrees*/
  void rotate( const PDBvector &, const PDBvector &, const double & );
                             /*rotate at new origin around axis some degrees*/
  void rotate( double **  );
  PDBatom rotateII( double ** );
  bool is_hydrogen( ) ;
  string get_resName_name( ) ;/*output residue name and atom name together*/
  string get_resName_nameII( ) ;/*correct for backbone atoms*/
  PDBvector normal_to_plane( PDBatom &, PDBatom & ) ; 
  double radius(PDBatom&,PDBatom&);

 protected:
  int serial ;         //PDBatom serial number
  string name ;           //PDBatom name
  string altLoc ;        //Alternate location indicator
  string resName ;    //Residue name
  string chainID ;        //Chain identifier
  int resSeq ;         //Residue sequence number
  string iCode ;         //Code for insertion of residues
  PDBvector r ;        //Orthogonal coordinates in Angstroms
  double occupancy ;   //Occupancy
  double tempFactor ;  //Temperature factor
  string segID ;   //Segment identifier, left-justified
  string element ; //Element symbol, right-justified
  string charge ;    //Charge on the PDBatom
  bool record ;     //true if coordinates are recorded

} ;


//----------- class nodePDBatom -----------


class nodePDBatom{
  
  friend class listPDBatom ;
  friend class PDBsalt ;

 public:
 nodePDBatom( const PDBatom & ) ;
 PDBatom & getTheAtom() ;
 nodePDBatom *getNext();

 protected:
  PDBatom theAtom ;
  nodePDBatom *next ;

} ;


//----------- class listPDBatom -----------


class listPDBatom{
  
  friend class PDBamino_acid ;
  friend ostream &operator<<(ostream & , const listPDBatom & ) ;/*declared
                                             also inside class nodePDBatom*/
  friend ifstream &operator>>(ifstream &, listPDBatom & ) ;/*declared also
                                                  inside class nodePDBatom*/

 public:
  listPDBatom( ) ;
  listPDBatom( const listPDBatom & ) ;
  listPDBatom( const char * ) ;/*contains the name of an input file*/
  ~listPDBatom( ) ;
  void empty( ) ;
  bool isEmpty( ) ;
  void insertAtFront( const PDBatom & );
  void rmFromFront( ) ;
  void insertAtBack( const PDBatom & );
  void insertAtBack( const listPDBatom & );
  int length( ) ;
  int importList( ifstream & ) ;
  int exportList(ostream &, int mode=0);
  PDBatom *getAtomAt( const int & ) ;
  PDBatom getAtomAtII( const int & ) ; /*same as *getAtomAt */
  string *printAtomAt( const int & ) ;
  bool matchResSeq( const PDBatom & ) ; 
  listPDBatom &operator=( listPDBatom &  ) ;
  listPDBatom &operator=( const listPDBatom &  ) ;
  PDBatom *heavyGeomCent( ) ;
  bool doContactHeavyAtom( PDBatom & , double & ) ;
  bool doContactHeavyAtom( PDBatom & , double & , int & ) ;
  bool doContactHeavyAtom( listPDBatom & , double & ) ;
  bool doContactAnyone( listPDBatom & , double & ) ;
  int numberAnyContacts( listPDBatom & , double & ) ;
  void changeResSeq( const int & ) ;
  void changeResSeqiCode( const int &, const string & ) ;
  void renumberFully( const int & ) ;
  int renumberFullyAtomSerialNumber( const int & ) ;
  void removeFullyiCodes();
  int removeRepeatedAtoms();
  bool isThereAtom( const string & ) ;
  PDBatom *getAtomFromName( const string & ) ;
  PDBatom *pointToAtomFromName( const string & ) ;
  void printChainID( string & ) ;
  void setChainID( string & ) ;
  void changeCoordAt( const int & , double *) ;
  void changeCoordAt( const int & , PDBvector & ) ;
  bool changeCoordOfAtomName( const string &, double * ) ;
  string * printTER( void ) ;
  double *getCoordOfAtomName( const string & ) ;
  int getSerialNumberOfAtomName( const string & ) ;
  char * printContiguousDist( void ) ;
  char * printNext_N_NeighDist( int & ) ;
  int createContactMap( ostream &, double & ) ;
  int createContactMap( double **, double & ) ;
  int createContactMap( listPDBatom &, double **, double & ) ;  
  double drms( listPDBatom &);
  double rmsd( listPDBatom &);
  double adapt_to_this_chain_by_rmsd( listPDBatom & ) ;
  double output_tranformed_chain_by_rmsd(listPDBatom &, listPDBatom & ) ;
  double rot_matrix_to_adapt_by_rmsd_to( listPDBatom &, double ** ) ;
  int diffContactMap( ifstream &, listPDBatom &, double &); 
  double NDrms( ifstream &,  listPDBatom & ) ;
  int numContacts( ifstream &, double & ) ;
  PDBvector get_CM( ) ;
  void translate( const PDBvector & ) ;
  void rotate( const PDBvector &, const double & ); /*rotate around axis some 
						      degrees*/
  void rotate( const PDBvector &, const PDBvector &, const double & );/*rotate at
                                              new origin around axis some degrees*/
  void rotate( double **  );
  double get_surrounding_sphere( ) ;
  double get_distance_between( int&, int& ) ;
  void get_coord_at_atom_number( int&, double* ) ;
  bool getCoordOfAtomWithAtomIndex( const int&, double * ) ;
  bool changeCoordOfAtomWithIndex( const int &, double * ) ;
  bool getNameOfAtomWithAtomIndex( const int &, string & ) ;
  PDBatom* returnAdressOfAtomWithAtomIndex( const int & ) ;
  double **dump_coordinates_to_array( ) ;
  int dump_coordinates_to_array2(double **) ;
  int dump_coordinates_to_chain( double **, int n=-1) ;
  void filterOnly( const string *, const int & ) ;
  double minD( listPDBatom & ) ;
  PDBatom getAtomWithAtomIndex( const int & ) ;
  bool hasAtomWithAtomIndex( const int & ) ;
  void remove_all_hydrogens( ) ;
  int extract(gen_atom_type, listPDBatom &);
  string connectRasmol();
  string red2blueRasmol();
  string blue2redRasmol();
  string white2blackRasmol();

 protected:
  nodePDBatom *firstPtr, *lastPtr ;
  nodePDBatom *getNewNode( const PDBatom & ) ;
 
} ;

//----------- class nodelistPDBatom -----------

class nodelistPDBatom{
  friend class listsPDBatom;

 public:
  nodelistPDBatom(const listPDBatom &, const string & = "");/*use default arg*/
  listPDBatom & getTheList() ;
  nodelistPDBatom *getNext();

 protected:
  string listID ;
  listPDBatom theList ;
  nodelistPDBatom *next ;
};

//------------ class listsPDBatom   --------------

class listsPDBatom{

  friend ostream &operator<<(ostream & , const listsPDBatom & ) ;/*declared
                                         also inside class nodelistPDBatom*/
  friend ifstream &operator>>(ifstream &, listsPDBatom & ) ;/*declared also
                                              inside class nodelistPDBatom*/

 public:
  listsPDBatom( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void empty( ) ;
  ~listsPDBatom( ) ;
  void insertAtBack( const listPDBatom &, const string & = "" );
  int length( ) ;
  listPDBatom &operator[]( const int & );  

 protected:
  nodelistPDBatom *firstPtr, *lastPtr ;
  nodelistPDBatom *getNewNode(const listPDBatom &, const string & = "");

};


//------------ class PDBamino_acid  --------------

class PDBamino_acid{

  friend ostream &operator<<(ostream & , const PDBamino_acid & ) ;
  friend ifstream &operator>>(ifstream &, PDBamino_acid &);
    

 public:
  PDBamino_acid( ) ;
  PDBamino_acid( string &) ;
  ~PDBamino_acid( ) ;
  bool isEmpty( ) ;
  void empty( ) ; 
  void assign( listPDBatom & , listPDBatom & ) ;
  void assign( const string & ) ;
  int exportAA(ostream &, int mode=0);
  string *printName( ) ;
  string getName( ) ;
  PDBamino_acid &operator=( PDBamino_acid & ) ;
  void addAtom( PDBatom & ) ;
  bool matchResSeq( PDBatom & ) ;
  string *printLastAtom( ) ;
  listPDBatom *getResidue( ) ;
  bool doContactResidues( PDBamino_acid & , double & ) ;
  bool doContactBackbones( PDBamino_acid & , double & ) ;
  bool doContactFirstResSecondBack( PDBamino_acid &, double & ) ;
  bool doContactFirstBackSecondRes( PDBamino_acid &, double & ) ;
  bool doContactResBack( PDBamino_acid & , double & ) ;
  bool doContactHN_HN( PDBamino_acid &, double &);
  bool doContactMETHYLC_METHYLC( PDBamino_acid &, double &);
  bool doContact( PDBamino_acid & , double & ) ;
  int  doContact( PDBamino_acid &, gen_atom_type, gen_atom_type, double &);
  int getResSeq( ) ;
  string getiCode();
  void changeFullyResSeq(const int &);
  void changeFullyResSeqiCode(int &, int &, string &);
  int renumberFullyAtomSerialNumber( const int &);
  void removeFullyiCodes();
  int removeRepeatedAtoms();
  bool is( const string & ) ;
  PDBatom *getAtomFromName( const string & ) ;
  PDBatom *pointToAtomFromName( const string &);
  PDBatom getAtomWithIndex( const int & ) ;
  bool hasAtomWithIndex( const int & ) ;
  void printChainID( string & ) ;
  bool hasAtomWithName( const string & ) ;
  bool changeCoordOfAtomName( const string &, double * ) ;
  bool changeCoordOfAtomWithIndex( const int &, double * ) ;
  double *getCoordOfAtomName( const string & ) ;
  void getCoordOfAtomName( const string &, PDBvector& );
  bool getCoordOfAtomWithAtomIndex( const int &, double * );
  bool getNameOfAtomWithAtomIndex( const int &, string & ) ;
  int numberOfAtoms( ) ;
  int getSerialNumberOfAtomName( const string & ) ;
  bool bonds_its_amide_H_to_the_backbone_O_of( PDBamino_acid & );
  void translate( const PDBvector & ) ;
  void rotate( const PDBvector &, const double & );/*rotate around axis some 
                                                                      degrees*/
  void rotate( const PDBvector &, const PDBvector &, const double & );/*rotate
				      at new origin around axis some degrees*/
  void rotate( double **  );/*rotation matrix*/
  PDBvector get_next_N_coords(void);
  PDBvector get_nextCA(void);
  void peptidePlaneAxis(PDBvector &, PDBvector &, PDBvector & );
  void intrinsicAxis( PDBvector &, PDBvector &, PDBvector & );
  bool align_intrinsic_axis_to_system_axis( ) ;
  void join_to( PDBamino_acid & ) ;
  void filterOnly( const string *, const int &n ) ;
  double minD(  PDBamino_acid & ) ;
  void remove_all_hydrogens( ) ;
  void fuseToSingleList( listPDBatom & ) ;
  int extract(gen_atom_type, listPDBatom &);
  PDBatom *generateHN(PDBatom *prevC) ;

 protected:

  string name ;
  int resSeq ;
  string iCode;
  listPDBatom backbone ;
  listPDBatom residue ;
} ;

bool (PDBamino_acid::*sel_doCont(cont_sel))(PDBamino_acid &,double &);/*return a pointer
  to appropriate doContact member function*/

class nodePDBamino_acid{

  friend class PDBchain ;

 public:
  nodePDBamino_acid( PDBamino_acid &);
  nodePDBamino_acid( PDBamino_acid & , const int &);
  ~nodePDBamino_acid();
  void changeFullyResSeq( const int & ) ;
  PDBamino_acid & getTheAA() ;
  nodePDBamino_acid *getNext();

 protected:
  int resSeq ;
  PDBamino_acid theAA ;
  nodePDBamino_acid *next ;
} ;



class PDBchain{

  friend ostream &operator<<(ostream & , const PDBchain & ) ;//declared
  //also inside class nodePDBoamino_acid
  friend ifstream &operator>>(ifstream &, PDBchain & ) ;//declared
  //also inside class nodePDBamino_acid

 public:
  PDBchain( ) ;
  ~PDBchain( ) ;
  PDBchain( ifstream & ) ;
  PDBchain( string & );
  PDBchain( char* );
  void empty( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void insertAtBack(PDBamino_acid &);
  void insertAtBack( PDBamino_acid & , const int & );
  void insertAtBack( PDBamino_acid & , const double, const double );
  int length( ) ;
  int importChain( ifstream & ) ;
  bool importChain(char *,const string &);
  int exportChain(char *, int mode=0, string endline="");
  PDBchain &operator=( PDBchain & ) ;
  int extract_segment(PDBchain &, int startResSeq, int endResSeq,
		      string startiCode=" ", string endiCode=" ");
  PDBamino_acid *getAminoFromIndex( const int & ) ;
  PDBamino_acid *pointToAminoFromIndex( const int & ) ;
  PDBamino_acid *getAminoFromResSeq( const int & ) ;
  nodePDBamino_acid *pointToNodeWithResSeqiCode(const int &,string ic=" ");
  string *printTER( ) ;
  string printChainID();
  int getIndexFromResSeq( const int & ) ;
  int getResSeqFromIndex( const int & ) ;
  string *printNameFromIndex( const int & ) ;
  string getNameFromIndex( const int & ) ;
  void renumberFully( const int & ) ;
  int renumberFullyAtomSerialNumber( const int & ) ;
  void removeFullyiCodes();
  int removeRepeatedAtoms();
  int removeResAtIndex(const int & );
  void simplify();
  bool isResNameAtIndex( const string &, const int & ) ;
  bool changeCoordOfAtomNameAtIndex( const string &, const int &, double * ) ;
  double *getCoordOfAtomNameAtIndex( const string &, const int & ) ;
  bool getCoordOfAtomNameAtIndex_II( const string &, const int &, double * ) ;
  bool getCoordOfAtomWithAtomIndex( const int &, double * ) ;
  int getSerialNumberOfAtomNameAtIndex( const string &, const int & ) ;/*Index
  refers to amino acid index along the chain*/
  PDBatom getAtomWithIndex( const int & ) ;
  bool changeCoordOfAtomWithIndex( const int &, double * ) ;
  bool changeCoordOfAtomWithIndex( const int &, PDBvector & ) ;
  PDBatom getAtomWithNameAndResIndex( const char *, const int & ) ;
  PDBatom *pointToAtomWithNameAndResIndex(const char *,const int &);
  bool getNameOfAtomWithAtomIndex( const int &, string & ) ;
  int numberOfAtoms( ) ;
  void createCAchain( listPDBatom & ) ;
  int createCAcontactMap( PDBchain &, double**, double & ) ;
  int createCBChain( listPDBatom & ) ;
  int createCBcontactMap( ostream &, double & ) ;
  int createCBcontactMap( double**, double & ) ;
  int createCBcontactMap( PDBchain &, double**, double & ) ;
  int createDistContMap( double**, double &, double & ) ;
  int createDistContMap( PDBchain &, double**, double &, double & ) ;
  int createHeavyAtomContactMap( int &, double**, double & ) ;
  int createSymmetricContactMap(int &, cont_sel, double **, double &);
  int createSymmetricContactMap(int &, gen_atom_type, double **, double &);
  int createAsymmetricContactMap(int &, cont_sel, double **, double &);
  int createAsymmetricContactMap(int &, gen_atom_type, gen_atom_type,
				 double **, double &);
  int createContactMap(int &, cont_sel, double **, double &);
  int createContactMap(int &, gen_atom_type,gen_atom_type,double **, double &);
  int createCACBchain( listPDBatom & ) ;
  double co(int &, cont_sel, double &);
  void insert_H( ) ;
  void remove_all_hydrogens( ) ;
  void translate( const PDBvector & ) ;
  void rotate( const PDBvector &, const double & );/*rotate around
						     axis some degrees*/
  void rotate( const PDBvector &, const PDBvector &, const double & );/*rotate
                                       at new origin around axis some degrees*/
  void rotate( double **  );
  void rotatePhiBy( const double &, const int & ) ;
  void rotatePsiBy( const double &, const int & ) ;
  double getPhiAtIndex( const int & ) ;
  double getPsiAtIndex( const int & ) ;
  double drmsPhiPsi( PDBchain & ) ;
  PDBvector get_CA_CM( ) ;
  PDBvector get_GC( ) ;/*GC stands for geometric center*/
  double get_CA_surrounding_sphere( ) ;
  void filterOnly( const string *, const int &n ) ;
  bool is_there_resName_name_at_resIndexII( const string &, const int & ) ;
  bool isAtomNameAtAminoAcidIndex(const string &, const int &);
  string output_one_letter_sequence_from_indexes( const int &, const int & ) ;
  string output_one_letter_sequence( ) ;
  void output_three_letter_sequence_from_indexes( const int &, const int &,
						  ostream & ) ;
  void output_three_letter_sequence( ostream & ) ;
  void fuseToSingleList( listPDBatom & ) ;
  double alignByRmsd( PDBchain &, const string *, const int & ) ;
  bool do_contact(PDBchain* ,bool (*do_aa_contact)(PDBamino_acid*,PDBamino_acid*));
  double average_graph_length(int &, cont_sel, double &);
  int number_of_contacts(int &, cont_sel, double &);
  void output_fasta_seq( ostream &, string & );
  protected:
  nodePDBamino_acid *firstPtr ;
  nodePDBamino_acid *lastPtr ;
  nodePDBamino_acid *getNewNode(PDBamino_acid &) ;
  nodePDBamino_acid *getNewNode( PDBamino_acid & , const int & ) ;
} ;




class nodePDBchain{

  friend class PDBchains ;

 public:
  nodePDBchain( PDBchain & , const string &) ;
  PDBchain & getTheChain() ;
  nodePDBchain *getNext();

 protected:
  string chainID ;
  PDBchain theChain ;
  nodePDBchain *next ;
} ;




class PDBchains{
  
  friend ostream &operator<<(ostream & , const PDBchains & ) ;//declared
  //also inside class nodePDBchain
  friend ifstream &operator>>(ifstream &, PDBchains & ) ;//declared
  //also inside class nodePDBchain

 public:
  PDBchains( ) ;
  PDBchains(char *);
  PDBchains(string &);
  ~PDBchains( ) ;
  void empty( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void insertAtBack( PDBchain & , const string & );
  int length( ) ;
  int importChains( ifstream & ) ;
  PDBchain* pointToPDBchainFromIndex( const int & );
  PDBchain* getPDBchainFromIndex( const int & );
  bool getPDBchainFromIndex(PDBchain &, const int & );
  PDBchain* pointToPDBchainFromChainID(string &);
  bool getPDBchainFromChainID(PDBchain &, const string & );

 protected:
  nodePDBchain *firstPtr ;
  nodePDBchain *lastPtr ;
  nodePDBchain *getNewNode( PDBchain & , const string & ) ;
} ;

/*===================================================================*/

#endif /*_PDBCLASSES_*/

