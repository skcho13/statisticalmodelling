* Encoding: windows-1252.
Matrix data variables=fam_back ability motivate courses achieve
    /contents=n corr sd mean
    /format=lower diagonal.



Begin data
1000 1000 1000 1000 1000
1.000
.417 1.000
.190 .205 1.000
.372 .498 .375 1.000
.417 .737 .255 .615 1.000
1.000 15.000 10.000 2.000 10.000
.000 100.000  50.000 4.000 50.000
End Data.

Regression  /Matrix=in(*)
        /variables= fam_back to achieve
        /descriptives=default,sig,n
        /stat=default ci
        /dep=achieve
            /method=enter fam_back to courses.

Regression  /Matrix=in(*)
        /variables= fam_back to achieve
        /stat=def ci
        /dep=courses
	/method=enter fam_back to motivate.

Regression  /Matrix=in(*)
        /variables= fam_back to achieve
        /stat=def ci 
        /dep=motivate
	/method=enter fam_back to ability.

Regression  /Matrix=in(*)
        /variables= fam_back to achieve
        /stat=def ci 
        /dep=ability
	/method=enter fam_back.
