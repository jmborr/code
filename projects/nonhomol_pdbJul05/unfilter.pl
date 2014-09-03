#!/usr/bin/perl

#pull only proteins for which Jeff calculated their homologous
chdir "/net/dell/01/users/skolnick/pdbhomoljul05";
@l=`ls -1 *.homol`;
foreach (@l){
    $id=substr($_,0,5);
    #find now length of the protein
    $_="/library/pdb_jul05/seq/$id.SEQ";
    if(-e $_){
	$_=`wc $_`;
	/(\d+)/; #get first number
	if($1<100){ printf "$id%4d\n", $1; }
    }
}

