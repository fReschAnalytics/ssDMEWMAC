$arlarray = @(2000)
$deltaarray = @(1.50,2.00,2.50,3.00)
$dimensionarray = @(2)
$ITERATIONS=10000
$lambdaarray = @(0.7)
$methodsarray = @("ssDMEWMC")
$waittimesarray = @(20,50,100,500)

for($d = 0; $d -lt $dimensionarray.length; $d++){
    for($m = 0; $m -lt $methodsarray.length; $m++){
        for($l = 0; $l -lt $lambdaarray.length; $l++){
            for($a = 0; $a -lt $arlarray.length; $a++){
                for($x = 0; $x -lt $deltaarray.length; $x++){
                    for($w = 0; $w -lt $waittimesarray.length; $w++){

                        $dimensionarray[$d]
                        $methodsarray[$m]
                        $lambdaarray[$l]
                        $arlarray[$a]
                        $deltaarray[$x]
                        $waittimesarray[$w]

                        docker run -it -v C:/Users/fresc/PycharmProjects/ssDMEWMAC/data:/out -e SIMULATION_TYPE="OOC" -e METH=$( $methodsarray[$m] ) -e DIME=$(  $dimensionarray[$d] ) -e ITERATIONS=10000 -e LAMBDA=$( $lambdaarray[$l] ) -e ARL=$( $arlarray[$a] ) -e DELTA=$( $deltaarray[$x] ) -e WAIT_TIME=$( $waittimesarray[$w] ) -e LOCAL_OR_GCP="LOCAL" ssdmewmac
                    }
                }
            }
        }
    }
}