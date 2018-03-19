function [ eta ] = tukey_window_create_old( width,alpha,y,dy )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


rm=width/2;

transition_region_width=rm-rm*(1-alpha);

y1=mean(y)-width/2-transition_region_width;
n1=round(y1/dy)+1;
n2=n1+round(transition_region_width/dy);
n3=n2+round(width/dy);
n4=n3+round(transition_region_width/dy);

eta=zeros(1,length(y));
r=y-mean(y);
eta(n1:n2)=.5-.5*cos(pi*(r(n1:n2)-rm*(1-alpha))./(alpha*rm) );
eta(n2:n3)=1;
eta(n3:n4)=.5-.5*cos(pi*(r(n3:n4)-rm*(1-alpha))./(alpha*rm) );

%eta(n3:n4)=fliplr(eta(n1:n2));

end

