--xdb.wrplab seq_demo seqdb_demo demo_pass;
SELECT prot.id, prot.seq, concat(sp.name," ",annot.descr)
 FROM prot,annot,sp
 WHERE prot.id=annot.prot_id and annot.db="sp" and annot.gi=sp.gi;
SELECT annot.prot_id, concat(sp.name," ",annot.descr)
 FROM annot,sp
 WHERE annot.prot_id=# AND annot.db="sp" and sp.gi=annot.gi;
SELECT prot.id,prot.seq  FROM prot where prot.id=#;
