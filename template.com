%chk={name}.chk
%nprocshared=4
%mem=4GB
# {method}/{basis} SCF=(VeryTight) Int=Ultrafine {opt}
Freq=VibRot Output=Pickett Geom=NoCrowd

{comment}

{charge} {multi}
{coords}

