#!/usr/bin/perl

`/bin/tar -zxvf nr.tar.gz`;
`/library/yzhang/blast/fastacmd -D >nr`;
`./pfilt nr >nr.filter`;
sleep(1);
`sync`;
#`/nfs/users/yzhang/bin/blast/formatdb -t $name\.filter -i $name\.filter`;
`/library/yzhang/blast/formatdb -i nr.filter -o T`;

exit();
