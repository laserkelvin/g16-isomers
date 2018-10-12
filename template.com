%chk={name}.chk
%nprocshared=4
%mem=4GB
# {method}/{basis} SCF=(QC,VeryTight) Int=Ultrafine {opt}
Freq=VibRot Output=Pickett

{comment}

{charge} {multi}
{coords}

