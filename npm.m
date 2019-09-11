function RHS = npm(t,R)

 global alpha beta mu gamma n0 max

RHS = zeros(size(R));
%%%%%%%%%%%%%%%%%%%%%
% Solve for moments%
%%%%%%%%%%%%%%%%%%%%%

dR(1) = alpha-mu*R(1)-2*beta*R(1)*R(2)+gamma*n0*(n0-1)*R(2); % x(1)=x(t) = the concentration of proteins in the normal conformation
dR(2) = -1*(mu+gamma*(2*n0-1))*R(2)+gamma*R(3); % x(2)=y(t)= the total number of aggregates
dR(3) = 2*beta*R(1)*R(2)-mu*R(3)-gamma*n0*(n0-1)*R(2); %x(3)=z(t) = the total number of prion protiens
% Note that aggregates are made of single prion protiens

RHS(1)=dR(1)';
RHS(2)=dR(2)';
RHS(3)=dR(3)';

%%%%%%%%%%%%%%%%%%%%%
% Solve for discrete%
%%%%%%%%%%%%%%%%%%%%%

RHS(4)=0;

%Soluble Protein
%RHS(1+3)= alpha - mu*R(1+3)-2*beta*R(1+3)*sum(R(9:max+3))+gamma*n0*(n0-1)*sum(R(9:max+3));

for i=n0:max
        %decay
        RHS(i)=RHS(i)-mu*R(i);
        
        %conversion
        RHS(i)= RHS(i)-2*beta*R(1)*R(i)+2*beta*R(1)*R(i-1);
       
        %fragmentation
        RHS(i)=RHS(i)-gamma*(i-1)*R(i);
        
%         for j= i+1:max
%             RHS(i)= RHS(i)+2*gamma*R(j);  
%         end
        for j=n0:i-1
            RHS(i)= RHS(i)-2*gamma*R(j);   
        end  
        RHS(i)= RHS(i)+2*gamma*R(2);
end
