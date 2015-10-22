function randnum = randfunc(xs,xe)

format long


rand('state', sum(1000*clock)); %%%Different random numer each time


  randnum=double(rand()) ;
  randnum=double((xs + ((xe-xs)*randnum)));

 	
	
