{

    if (step%natoms == 0) 

    {

        print natoms
        print int(step/natoms) + 1
        print "H", $1*acell*0.529177, $2*acell*0.529177, $3*acell*0.529177

    }

    else print "H", $1*acell*0.529177, $2*acell*0.529177, $3*acell*0.529177
    step+=1

}
