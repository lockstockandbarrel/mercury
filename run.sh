#!/bin/bash
(
exec 2>&1
fpm @clean
fpm @run
) | tee $(basename $0 .sh).log
echo 'Summary'
grep 'runs for' $(basename $0 .sh).log
exit

wget "ftp://ftp.urbanjost.altervista.org/REMOVE/timing1.tgz"

