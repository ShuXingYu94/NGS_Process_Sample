Tmp_fifo=./$$.fifo
mkfifo $Tmp_fifo
exec 6<> $Tmp_fifo
rm -rf $Tmp_fifo

for i in `seq 50`     # Threds number
do

  echo >&6

done

for i in {1..254}
do

  read -u 6
  {

    echo "$i"         # Main procedure
    sleep 1

  echo >&6
  }&
done
wait
exec 6>&-
echo 'done'