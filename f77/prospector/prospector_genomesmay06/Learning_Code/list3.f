c       divide all queries into packs of ten. For each pack, create a
c       LIST.targ file containing the query identifiers and a nametarg
c       file containing a prefix of the form pdbbX or pdbbXX or pdbXXX

	character*255 namep(25000),name
	character*2 namefiles(99)
	character*3 namefiles3(300)	
	character*4 pdbf
	character*3 pdb
	character*1 namefiles1(9)
        data namefiles1/'1','2','3','4','5','6','7','8','9'/
	data namefiles/'_1','_2','_3','_4','_5','_6','_7','_8','_9',
     &	'10','11','12','13','14','15','16','17','18','19',
     &	'20','21','22','23','24','25','26','27','28','29',
     &	'30','31','32','33','34','35','36','37','38','39',
     &	'40','41','42','43','44','45','46','47','48','49',
     &	'50','51','52','53','54','55','56','57','58','59',
     &	'60','61','62','63','64','65','66','67','68','69',
     &	'70','71','72','73','74','75','76','77','78','79',
     &	'80','81','82','83','84','85','86','87','88','89',
     &	'90','91','92','93','94','95','96','97','98','99'/
	data namefiles3/'100',
     &  '101','102','103','104','105','106','107','108','109',
     &	'110','111','112','113','114','115','116','117','118','119',
     &	'120','121','122','123','124','125','126','127','128','129',
     &	'130','131','132','133','134','135','136','137','138','139',
     &	'140','141','142','143','144','145','146','147','148','149',
     &	'150','151','152','153','154','155','156','157','158','159',
     &	'160','161','162','163','164','165','166','167','168','169',
     &	'170','171','172','173','174','175','176','177','178','179',
     &	'180','181','182','183','184','185','186','187','188','189',
     &	'190','191','192','193','194','195','196','197','198','199',
     &  '200',
     &  '201','202','203','204','205','206','207','208','209',
     &	'210','211','212','213','214','215','216','217','218','219',
     &	'220','221','222','223','224','225','226','227','228','229',
     &	'230','231','232','233','234','235','236','237','238','239',
     &	'240','241','242','243','244','245','246','247','248','249',
     &	'250','251','252','253','254','255','256','257','258','259',
     &	'260','261','262','263','264','265','266','267','268','269',
     &	'270','271','272','273','274','275','276','277','278','279',
     &	'280','281','282','283','284','285','286','287','288','289',
     &	'290','291','292','293','294','295','296','297','298','299',
     &  '300',
     &  '301','302','303','304','305','306','307','308','309',
     &	'310','311','312','313','314','315','316','317','318','319',
     &	'320','321','322','323','324','325','326','327','328','329',
     &	'330','331','332','333','334','335','336','337','338','339',
     &	'340','341','342','343','344','345','346','347','348','349',
     &	'350','351','352','353','354','355','356','357','358','359',
     &	'360','361','362','363','364','365','366','367','368','369',
     &	'370','371','372','373','374','375','376','377','378','379',
     &	'380','381','382','383','384','385','386','387','388','389',
     &	'390','391','392','393','394','395','396','397','398','399'/

     	open(unit=1,file='input') !input contains just a word to be used as prefix
	read(unit=1,fmt='a4')pdbf
     	open(unit=2,file='input2') !input2 contains just a word 
	read(unit=2,fmt='a3')pdb
	open(unit=6,file='out')
	open(unit=11,file='namesize') !namesize contains length of word identifying each template
	read(11,*)nsize
	close(11)
	open(unit=10,file='LIST') !list of queries
	ics=0
	read(10,*)naln !number of queries
	do id=1,naln
	   read(10,1)name(1:nsize) !read one query identifier at a time
	   ics=ics+1
	   namep(ics)=name(1:nsize) !store names of queries on namep
	end do

	if(mod(naln,10).ne.0)then
	   ncs=naln/10 +1           !divide queries into packs of ten
	else
	   ncs=naln/10
	end if
	write(6,*)ncs
	open(unit=7,file='nruns')
	write(7,*)ncs       
	do imv=1,ncs !create files for each pack of then queries
	   if(imv .gt.9 .and. imv .lt.100)then
	      OPEN(unit=30,file='nametarg'//namefiles(imv)) !nametargXX
	   elseif(imv .lt.10)then
	      OPEN(unit=30,file='nametarg'//namefiles1(imv)) !nametargX
	   end if
	   if(imv .lt. 100)then
	      OPEN(unit=20,file='LIST.targ'//pdbf//namefiles(imv)) !LIST.targpdbbXX
	      write(30,22)pdbf,namefiles(imv) !the nametarg file contains "pdbfXX"
 22	      format(a4,a2)
	   else
	      iimv=imv-99	
	      OPEN(unit=30,file='nametarg'//namefiles3(iimv)) !nametargXXX
	      OPEN(unit=20,file='LIST.targ'//pdb//namefiles3(iimv)) !LIST.targpdbXXX
	      write(30,210)pdb,namefiles3(iimv) !the nametarg file contains pdbXXX, not pdbbXXX
 210	      format(a3,a3)
	   end if
	   if(imv.lt.ncs)then !number of queries in the pack. Last pack may have less than ten
	      id=10
	   else
	      id=naln-(ncs-1)*10
	   end if
	   write(20,2)id,pdbf !first line in the LIST.targ file contains "10 pdbb"
 2	   format(i3,1x,a4)
	   do ii=1,id
	      write(20,1)namep(ii+10*(imv-1)) !subsequent lines contain query identifiers
	   end do
 1	   format(a)
	end do
	close(20)
	stop
	end 