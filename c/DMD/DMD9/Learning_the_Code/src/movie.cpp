#include "movie.h"

/* read the header of the movie file */
int movie::read_header(){
  char * password ="BCPkino";
  char * password1="BCPkina";
  l_pass=strlen(password);
  
  int i, j;
  
  /*passwd and version comparison*/
  char * buf = (char *)malloc(l_pass+1);
  fread(buf, sizeof(char), l_pass, mov_path);
  buf[l_pass]='\0';
  version=0;
  if(strcmp(password,buf)){
    version=1;
    if(strcmp(password1,buf)){
      errMsg = "wrong movie format.";
      return 0;
    }
  }
  free(buf);
  
  /*movie frames*/
  if(fread_bo(&nFrames, sizeof(nFrames), 1, mov_path) < 1) 
    return 0;
  if(DEBUG)printf("number of frames=%ld\n",nFrames);
  frameOffset.push_back(0);
  new_type.push_back(1);
  new_bond.push_back(1);  
  
  /*system size components*/
  if(fread_bo(box, sizeof(double), 3, mov_path) < 3)
    return 0;
  if(DEBUG)printf("X=%lf Y=%lf Z=%lf\n",box[0],box[1],box[2]);
  
  /* delta T */
  if(fread_bo(&dt, sizeof(double), 1, mov_path) < 1) 
    return 0;
  if(DEBUG)printf("dt=%lf\n",dt);
  
  /* number of diff. types */
  if(fread_bo(&nAtomTypes, sizeof(nAtomTypes), 1, mov_path) < 1) 
    return 0;
  if(DEBUG)printf("nat=%d\n",nAtomTypes);
  
  atomColors=(short *)malloc(nAtomTypes*sizeof(short));
  atomRadii =(double *)malloc(nAtomTypes*sizeof(double));
  /* types of atoms */
  if(fread_bo(atomColors, sizeof(short), nAtomTypes, mov_path) < nAtomTypes) 
    return 0;
  /* radius */
  if(fread_bo(atomRadii, sizeof(double), nAtomTypes, mov_path) < nAtomTypes)
    return 0;
  
  /*number of atoms*/
  if(fread_bo(&nAtoms, sizeof(nAtoms), 1, mov_path) < 1) 
    return 0;
  if(DEBUG)printf("n_atom=%d\n",nAtoms);
  
  /*corrdinate buffer*/
  r=(double **)malloc(sizeof(double *)*nAtoms);
  r[0]=(double *)malloc(sizeof(double)*nAtoms*3);
  for(i=1;i<nAtoms;i++)
    r[i]=r[i-1]+3;
  
  /*list of atoms types */
  atomTypes=(short *)malloc(nAtoms*sizeof(short));
  if(fread_bo(atomTypes, sizeof(short), nAtoms, mov_path) < nAtoms) 
    return 0;  
  
  /* list of bonds*/
  /*---------------FORMAT----------------
    0, ..., 0 (... means the index of friends that are beger than i)
    1, ..., 1
    ...
    ...
    ---------------END-------------------*/
  bonds = (short**)malloc(nAtoms*sizeof(short*));
  bonds[0] = (short*)malloc((nAtoms*(nAtoms+1))/2*sizeof(short));
  for(i=1; i<nAtoms; i++)
    bonds[i]=bonds[i-1]+(nAtoms+1 - i );
  if(!read_bond()) return 0;
  
  /* number of params */ 
  if(fread_bo(&nParameters, sizeof(nParameters), 1, mov_path) < 1) 
    return 0;
  if(DEBUG)printf("n_param=%d\n",nParameters);
  parameterValues=(double*)malloc(nParameters*sizeof(double));
  parameterNames = (char**)malloc(nParameters*sizeof(char*)); 
  
  /*list of Parameters*/
  for(i=0; i<nParameters; i++){
    if(fread_bo(&s_dummy, sizeof(s_dummy), 1, mov_path) < 1) 
      return 0;
    parameterNames[i] = (char*)malloc((s_dummy+1));
    if(fread(parameterNames[i], s_dummy, 1, mov_path) < 1)
      return 0;
    parameterNames[i][s_dummy]='\0';
    if(DEBUG)printf("%s\n",parameterNames[i]);
  }
  
  if(box[2]!=0)nDimensions = 3;
  else nDimensions = 2;
  
  for(i=0; i<nDimensions; i++)
    scaleFactors[i]=box[i]/(double)65536;
  
  return ftell(mov_path);
}

/*in the case of the corrupt binary file,
  create an array of the offsets for each valid frame, so that
  the search of different frames becomes easy.*/
void movie::make_offset(){
  long begin=ftell(mov_path);
  int i;
  int good=1;
  long offset;
  while(good){
    offset = ftell(mov_path);
    /*new list of atoms types ?*/
    if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) break;
    /*yes--------read the new type list*/
    if(c_dummy){
      new_type.push_back(1);
      if(fseek(mov_path, nAtoms*sizeof(short), SEEK_CUR)) break;
    }
    else new_type.push_back(0);
    
    /*new list of bonds ?*/
    if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1 ) break;
    /*yes--------read the bond list*/
    if(c_dummy){
      new_bond.push_back(1);
      for(i=0; i<nAtoms; i++){
	if(fread_bo(&s_dummy, sizeof(s_dummy), 1, mov_path) < 1) {
	    good=0;
	    break;
	}
	if(s_dummy!=i){
	  errMsg = "wrong format of bonds";
	  good=0;
	  break;
	}
	while(1){
	  if(fread_bo(&s_dummy, sizeof(s_dummy), 1, mov_path) < 1){
	    good=0;
	    break;
	  }
	  if(s_dummy==i) break;
	}
	if(!good) break;
      }
      if(!good) break;
    }
    else new_bond.push_back(0);
    
    /* coords + params + EOFrame*/
    if(fseek(mov_path, nAtoms*nDimensions*sizeof(short)
	     + nParameters*sizeof(double) + sizeof(short),
	     SEEK_CUR)) break;
    
    readed++;
    frameOffset.push_back(offset);
    if(feof(mov_path)) break;
  }
  /*go back to the previous point*/
  fseek(mov_path, begin, SEEK_SET);
  
  if(DEBUG)printf("Fames: %ld %ld %ld\n",frameOffset.size()-1, readed ,nFrames);
  if(readed!=nFrames){
    errMsg += " incomplete movie file.";
    nFrames = readed;
  }
}

int movie::readChanges(int type_frame, int bond_frame){
  int bond_readed = 0;
  if(type_frame >= 0){
    fseek(mov_path, frameOffset[type_frame], SEEK_SET);
    if(type_frame==0){/*read from the header*/
      int skip;
      skip = l_pass + sizeof(nFrames) + 3*sizeof(double) + sizeof(double)
	+ sizeof(nAtomTypes) + nAtomTypes*(sizeof(short) + sizeof(double))
	+ sizeof(nAtoms);
      fseek(mov_path, skip, SEEK_CUR);
      if(fread_bo(atomTypes, sizeof(short), nAtoms, mov_path) < nAtoms) return 0;
      if(bond_frame == type_frame){
	if(!read_bond()) return 0;
	bond_readed = 1;
      }
    }
    else{
      /*new list of atoms types ?*/
      if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
      if(c_dummy){/*read the new type list*/
	if(fread_bo(atomTypes, sizeof(short), nAtoms, mov_path)< nAtoms) return 0; 
      }
      else {errMsg += "no newTypes!"; return 0;}
      if(bond_frame==type_frame){
	if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
	if(c_dummy){
	  if(!read_bond()) return 0;
	}
	else {errMsg += "no newBonds!"; return 0;}
	bond_readed = 1;
      }
    }
  }
  
  if(bond_frame >=0 && !bond_readed){
    fseek(mov_path, frameOffset[bond_frame], SEEK_SET);      
    if(bond_frame==0){
      int skip;
      skip = l_pass+sizeof(nFrames)+3*sizeof(double)+sizeof(double)
	+sizeof(nAtomTypes)+nAtomTypes*(sizeof(short)+sizeof(double))
	+sizeof(nAtoms)+nAtoms*sizeof(short);
      fseek(mov_path, skip, SEEK_CUR);
      if(!read_bond()) return 0;
    }
    else{
      /*new list of atoms types ?*/
      if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path)<1) return 0;
      if(c_dummy)/*read the new type list*/
	fseek(mov_path, nAtoms*sizeof(short), SEEK_CUR);
      
      if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
      if(c_dummy){
	if(!read_bond()) return 0;
      }
      else{
	errMsg += "no newBonds!"; 
	return 0;
      }
    }
  }
  return 1;
}

/*read current frame started by current file poonter*/
int movie::read_frame(char& typeChanges, char& bondChanges){
  int i,j;
  unsigned short coord;
  
  /*new list of atoms types ?*/
  if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path)<1) return 0;
  typeChanges = c_dummy;
  if(c_dummy)/*read the new type list*/
    if(fread_bo(atomTypes, sizeof(short), nAtoms, mov_path)< nAtoms) return 0;  
  
  /*new list of bonds ?*/
  if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
  bondChanges = c_dummy;
  if(c_dummy) if(!read_bond()) return 0;
  
  /* coords */
  if(version)
    for(i=0; i<nAtoms; i++)
      for(j=0; j<nDimensions; j++){
	if(fread_bo(&coord, sizeof(coord), 1, mov_path) < 1) return 0;
	r[i][j]=(coord)*scaleFactors[j];
      }
  else
    for(i=0; i<nAtoms; i++)
      for(j=0; j<nDimensions; j++){
	if(fread_bo(&coord, sizeof(coord), 1, mov_path)<1) return 0;
	r[i][j]=coord;
      }
  
  /*params */
  if(fread_bo(parameterValues, sizeof(double), nParameters, mov_path) < nParameters) 
    return 0;
  if(DEBUG)
    for(i=0; i<nParameters; i++)
      printf("%s: %20.18lf\n", parameterNames[i], parameterValues[i]);
  
  /*END of FRAME*/
  do{
    if(fread_bo(&s_dummy,sizeof(s_dummy), 1 ,mov_path)<1) return 0;
  }while(s_dummy!=0);
  
  return 1;
}

movie::movie(char* movie_name, int DEBUG_OPT){
  /*init*/
  errStatus = 0;
  mov_path = NULL;
  nFrames = 0;
  nParameters =0;
  parameterNames = NULL;
  parameterValues = NULL;
  atomColors = NULL;
  atomRadii = NULL;
  atomTypes = NULL;
  r = NULL;
  
  if(movie_name){
    //cout << movie_name << endl;
    mov_path = fopen(movie_name,"rb");
  }
  if(!mov_path){
    errStatus = 1;
    errMsg = "Can not open the binary movie file " + string(movie_name);
    return;
  }
  DEBUG = DEBUG_OPT;
  bo = byte_order();
  if( !(read_header()) ){
    errStatus = 2;
    errMsg = "Can not read the movie header: " + errMsg;
    fclose (mov_path);
    return;
  }
  readed=0;
  if(!nFrames)
    make_offset();
  else if(nFrames<0)
    make_offset();
  current = 0;
}

movie::~movie(){
  if(mov_path){
    fclose(mov_path);
    free(atomColors);
    free(atomRadii);
    free(atomTypes);
    free(parameterValues);
    frameOffset.clear();
    new_type.clear();
    new_bond.clear();
    for(int i=0; i<nParameters; i++)
      free(parameterNames[i]);
    free(parameterNames);
    if(r)free(r[0]);
    free(r);
  }
}

int movie::nextFrame(){
  char newType, newBond;
  long top = 0;
  if( (readed<nFrames) && current == readed) top = ftell(mov_path);
  if(read_frame(newType, newBond)){
    if(top){
      readed ++ ;
      frameOffset.push_back(top);
      new_type.push_back(newType);
      new_bond.push_back(newBond);
    }
    current++;
    return 0;
  }
  else{
    return 1;
  }
}

int movie::prevFrame(){
  char newType, newBond;
  if(current>1){
    if((!new_type[current] && !new_bond[current])||
       new_type[current-1] && new_bond[current-1])
      {
	fseek(mov_path, frameOffset[current-1], SEEK_SET);
	read_frame(newType, newBond);
	current--;
      }
    else{
      int type_frame, bond_frame;
      for(int i=current-1; i>=0; i--){
	if(new_type[i]){
	  type_frame=i;
	  break;
	}
      }
      for(int i=current-1; i>=0; i--){
	if(new_bond[i]){
	  bond_frame=i;
	  break;
	}
      }
      if(!readChanges(type_frame, bond_frame)) return 1;
      fseek(mov_path, frameOffset[current-1], SEEK_SET);
      read_frame(newType, newBond);
      current--;
    }
    return 0;
  }
  else{
    return 1;
  }
}

int movie::jumpto(int target){
  if(target <= nFrames && target >0){
    if(target <= readed){
      if(!new_type[target] || !new_bond[target]){
	/*do we need to reread the bond and types lists?*/
	int type_frame=-1;
	int bond_frame=-1;
	if(!new_type[target]){
	  if(target>current){
	    for(int i=current+1; i<target; i++){
	      if(new_type[i]){
		type_frame = i; break;
	      }
	    }
	    if(type_frame>0){
	      for(int i=target-1; i>current; i--){
		if(new_type[i]){
		  type_frame= i; break;
		}
	      }
	    }
	  }
	  else if(target<current){
	    for(int i=target+1; i<current; i++){
	      if(new_type[i]){
		type_frame = i; break;
	      }
	    }
	    if(type_frame>0){
	      for(int i=target-1; i>=0; i--){
		if(new_type[i]){
		  type_frame = i; break;
		}
	      }
	    }
	  }
	}
	if(!new_bond[target]){
	  if(target>current){
	    for(int i=current+1; i<target; i++){
	      if(new_bond[i]){
		bond_frame = i; break;
	      }
	    }
	    if(bond_frame>0){
	      for(int i=target-1; i>current; i--){
		if(new_bond[i]){
		  bond_frame = i; break;
		}
	      }
	    }
	  }
	  else if(target<current){
	    for(int i=target+1; i<current; i++){
	      if(new_bond[i]){
		bond_frame = i; break;
	      }
	    }
	    if(bond_frame>0){
	      for(int i=target-1; i>=0; i--){
		if(new_bond[i]){
		  bond_frame = i; break;
		}
	      }
	    }
	  }
	}
	
	if(!readChanges(type_frame, bond_frame)) return 0;
      }
      if(target!=current){
	fseek(mov_path, frameOffset[target], SEEK_SET);
	char c_dummy1, c_dummy2;
	read_frame(c_dummy1, c_dummy2);
	current = target;
      }
    }
    else{
      long offset;
      for(int i=current+1; i<target; i++){
	offset = ftell(mov_path);
	/*new list of atoms types ?*/
	if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
	/*yes--------read the new type list*/
	if(c_dummy){
	  new_type.push_back(1);
	  if(fread_bo(atomTypes, sizeof(short), nAtoms, mov_path)< nAtoms) return 0;  
	}
	else new_type.push_back(0);
	
	/*new list of bonds ?*/
	if(fread_bo(&c_dummy, sizeof(c_dummy), 1, mov_path) < 1) return 0;
	if(c_dummy){
	  new_bond.push_back(1);
	  if(!read_bond()) return 0;
	}
	else new_bond.push_back(0);
	
	/* coords + params + EOFrame*/
	if(fseek(mov_path, nAtoms*nDimensions*sizeof(short)
		 + nParameters*sizeof(double) + sizeof(short),
		 SEEK_CUR)) return 0;
	
	readed++;
	current++;
	frameOffset.push_back(offset);
      }
      char c_dummy1, c_dummy2;
      offset = ftell(mov_path);
      if(!read_frame(c_dummy1, c_dummy2))return 0;
      new_type.push_back(c_dummy1);
      new_bond.push_back(c_dummy2);
      frameOffset.push_back(offset);
      readed++;
      current++;
    }
    return 1;
  }
  else return 0;
}


/*
int main(int argc, char** argv){
  if(argv[1]){
    movie m(argv[1],1);
    if(m.getErrorStatus()){
      cout << m.getMsg() << endl;
    }
    int i;
    while(1){
      cout << " Total: " << m.getNFrames() << endl;
      cout << " readed:" << m.getReaded() << endl;
      cout << " size: " << m.getSize() << endl;
      cout << " Current:" << m.getCurrentFrame() << endl;
      cout << " input : " << endl;
      cout << " 1. goon; 2. jump; " << endl;
      
      cin >> i;
      switch(i){
      case 1:
	m.nextFrame();
	break;
      case 2:
	int t;
	cout << "Where to jump" << endl;
	cin >> t;
	m.jumpto(t);
	break;
      }
    }
  }
}

*/
