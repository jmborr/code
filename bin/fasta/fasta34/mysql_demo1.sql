localhost database user password;
SELECT proteins.gid, proteins.sequence FROM proteins,gi2taxid, nested WHERE proteins.gid = gi2taxid.gid and gi2taxid.taxid= nested.taxid and nested.leftid between 198215 and 203126 ;
select proteins.gid, concat(proteins.gid," ",proteins.description) from proteins where proteins.gid=#;
select gid, sequence from proteins where gid=#;
