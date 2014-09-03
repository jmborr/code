#!/usr/bin/perl

$db='/net/dell/01/users/skolnick/pdbhomoljul05';
open(in,"/net/ibm/03/jose/Projects/nonhomol_pdbJul05/run/unfiltered.list");
$l=''; while(<in>){ chomp; $l .= $_; } close(in);
$id=substr($l,0,5,''); #pull and remove $id from the front of $l
$done='';
while(!($done=~$id)){
    $done.=$id; #mark $id as dealt with
    $l.=$id; #add $id to the end of $l
    open(in,"$db/$id.homol") or die "can't find $id.homol  \n";
    while(<in>){
	$id2=substr($_,0,5);
	$l=~s/$id2//; #remove homologous from list
    }
#    print length($l)."\n";
    close(in);
    $id=substr($l,0,5,''); #pull and remove next entry
#    print "id=$id done=$done\n";
}
while(length($done)>0){
    print substr($done,0,5,'')."\n";
}
