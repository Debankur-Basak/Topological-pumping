function output=AAH_Model_Hamiltonian(N,Jab,Jbc,Jca,delta,pbc)

NN=3*N;                              % Each cell has three sites
output=zeros(NN);

beta_A=delta;                        % Onsite energy for site a
beta_B=delta;                        % Onsite energy for site b
beta_C=delta;                        % Onsite energy for site c


%% Diagonal elements of the matrix.

for ii=1:3:NN
    output(ii,ii)=beta_A;
    output(ii+1,ii+1)=beta_B;
    output(ii+2,ii+2)=beta_C;

end

%% Off-diagonal elements of the matrix.
for ii=1:3:NN

    output(ii,ii+1)=Jab;
    output(ii+1,ii)=conj(Jab);
    
    if ii<NN

        output(ii+1,ii+2)=Jbc;
        output(ii+2,ii+1)=conj(Jbc);

    end

    if ii<3*(N-1)

        output(ii+2,ii+3)=Jca;
        output(ii+3,ii+2)=conj(Jca);

    end
end

%% Existance of Periodic Boundary Condition.

if pbc==1
    output(1,NN)=conj(Jca);
    output(NN,1)=Jca;
end

end