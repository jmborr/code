localhost database user password;
SELECT proteins.gid, proteins.sequence FROM proteins,swissprot WHERE proteins.gid=swissprot.gid;
select proteins.gid, concat(swissprot.spid," ",proteins.description) from proteins,swissprot where proteins.gid=# AND swissprot.gid=proteins.gid;
select gid, sequence from proteins where gid=#;
