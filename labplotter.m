function resultvector=labplotter(M,N,Nchan,lim)
	serv=5*M/(N*Nchan);
	l=3*serv;
	ro=l/serv;
	P0a=1;
	for q=1:Nchan
		P0a=P0a*ro/q;
	endfor
	P0a=P0a*ro/(Nchan-ro);
	P0base=1;
	P0b=1;
	for q=1:Nchan
		P0base=P0base*ro/q;
		P0b=P0b+P0base;
	endfor
	P0=1/(P0a+P0b);
	k=1:lim;
	P=arrayfun(@(x) P0*(ro^x)/(gamma(Nchan+1)*Nchan^(x-Nchan)),k);
	busy=P0*(ro^Nchan)/(gamma(Nchan)*(Nchan-ro));
	queued=busy*ro/N;
	base=1;
	value=1;
	for q=1:Nchan-1
		base=base*ro/q;
		value=value+base;
	endfor
	meanqueue=P0*(ro*value+(ro^(Nchan+1))*(Nchan+1-ro)/(gamma(Nchan)*((Nchan-ro)^2)));
	meanlength=(ro^(Nchan+1)*P0)/(gamma(Nchan)*(Nchan-ro)^2);
	basic=1/gamma(Nchan+1);
	meanfree=0;
	for q=1:Nchan
		basic=basic*ro*(Nchan-q+1);
		meanfree=meanfree+q*basic;
	endfor
	meanfree=meanfree*P0;
	meanbusy=Nchan-meanfree;
	meanwait=P0*(ro^Nchan)/(serv*gamma(Nchan)*(Nchan-ro)^2);
	wait=meanwait*ro*serv;
	meantime=meanwait+1/serv;
	sumtime=wait+ro;
	plot(k,P);
	resultvector=[queued busy meanqueue meanlength meanfree meanbusy meanwait wait meantime sumtime];
endfunction
