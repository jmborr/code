#!/usr/bin/perl -w
###########################################################
#
#DISpro: Protein Disorder region Prediction Programs
#configure.pl: to configure the installation of DISpro
#
#Author: Jianlin cheng
#Date: July 4, 2004
#Institute for Genomics and Bioinformatics
#School of Information and Computer Science
#University of California, Irvine
#
##########################################################

#######Customize settings here##############
#
#set installation directory of DISpro1.0
$install_dir = "/gpfs1/active/jose/code/bin/dispro1.0/";

#set the fullpath of sspro4 installation directory
$pspro_dir = "/gpfs1/active/jose/code/bin/sspro4/";
#
######End of Customization##################


################Don't Change the code below##############
if (! -d $install_dir)
{
	die "can't find installation directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}
if (! -d $pspro_dir)
{
	die "can't find Pspro directory.\n";
}
if ( substr($pspro_dir, length($pspro_dir) - 1, 1) ne "/" )
{
	$pspro_dir .= "/"; 
}

#check if the installation directory is right
#the configuration file must run in the installation directory
$cur_dir = `pwd`;  
chomp $cur_dir; 
$configure_file = "$cur_dir/configure.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
	die "Please check the installation directory setting and run the configure program there.\n";
}

$bin_dir = "${install_dir}bin/";
$model_dir = "${install_dir}model/";
$script_dir = "${install_dir}script/";
$server_dir = "${install_dir}server/";
$test_dir = "${install_dir}test/";

if ( ! -d $bin_dir || ! -d $model_dir || ! -d $script_dir
   || ! -d $server_dir || ! -d $test_dir )
{
	die "some sub directories don't exist. check the installation tar ball.\n";
}

#generate model definition file
$diso_exe = "${server_dir}predict_seq";
$diso_sh = "${bin_dir}predict_diso.sh";
$diso_model_dir = "${model_dir}diso/";
$diso_model_def = "${model_dir}model.def";
$pspro_exe = "${pspro_dir}bin/predict_ss_sa.sh"; 

print "generate diso server script...\n";
open(SERVER_SH, ">$diso_sh") || die "can't write diso shell script.\n";
print SERVER_SH "#!/bin/sh\n#predict disorder region for one sequence.\n";
print SERVER_SH "if [ \$# -ne 2 ]\n";
print SERVER_SH "then\n\techo \"need 2 parameters:fasta seq_file, output prefix.\"\n\texit 1\nfi\n";
print SERVER_SH "${script_dir}predict_diso.pl $pspro_exe $diso_exe $diso_model_def \$1 \$2 \n"; 
close SERVER_SH;

print "generate diso model definition file...\n";
opendir(MODEL_DIR, "$diso_model_dir") || die "can't open the diso model dir.\n";
open(MODEL_DEF, ">$diso_model_def") || die "can't create diso model def file.\n";
@file_list = readdir(MODEL_DIR);
closedir(MODEL_DIR); 
if (@file_list < 3)
{
	die "can't find diso model file. check intallation tar ball.\n";
}
$model_num = @file_list;
$model_num -= 2; 
print MODEL_DEF "$model_num 2\n";
while (@file_list)
{
	$model_file = shift @file_list;
	if ($model_file ne "." && $model_file ne ".." && $model_file ne "CVS")
	{
		print MODEL_DEF "$diso_model_dir$model_file\n";
	}
}
close MODEL_DEF; 
`chmod 755 $diso_sh`; 


