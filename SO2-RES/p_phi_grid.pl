#!/usr/bin/perl

$istate = $ARGV[0];
$axes = $ARGV[1];
$sigma = $ARGV[2];

$j    = 40;
$dp   = 1.0;
$dphi = 1.0;

$np=0;
#foreach $sigma (0..1)
{
  for ($p=-$j; $p<=$j; $p=$p+$dp)
  {
    for ($phi=0.0; $phi<=180.0; $phi=$phi+$dphi)
    {
      ++$np;

      #print "$np --- run $sigma $p $phi\n";
      #`./SO2_ames2_p12.run $j $sigma $p $phi`;

      if ($axes==1)
      {
        $fname = "so2_pcossin_j".$j."_s".$sigma."_p".$p."_phi".$phi.".out";
      }
      else
      {
        $fname = "so2_sincosp_j".$j."_s".$sigma."_p".$p."_phi".$phi.".out";
      }

      $enr = `./read_enr.pl $fname $istate`;
      printf("  %12.6f  %12.6f  %12.6f\n", $p, $phi, $enr);
    }
  }
}
#print "$np\n";
