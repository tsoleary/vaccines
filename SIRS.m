function dSIR = SIRS(t,SIR, alpha, beta, gamma, delta, epsilon, contact, immunity)
    SIR=SIR(1:end);
    n=size(contact,1);
    N=sum(SIR);
    S=SIR(1);
    I=SIR(2:n+1)';
    R=SIR(n+2:end)';
    immunity=ones(size(immunity))-immunity;
    %initialize rates
    dSdt = 0;
    dIdt(1:n) = 0;
    dRdt(1:n) = 0;
    DELTA(1:n,1:n)=0;
    NBOR(1:n,1:n)=0;
    for i=1:n
        for j=1:n
            if delta==0 %disabled transcending immunity case
                
                %allow immunity only for current strain   
                if contact(i,j)==0 %prior error: all were set to 0, instead of just distance=0
                    DELTA(i,j)=0;
                else
                    DELTA(i,j)=1; %full infection rate for others
                end
                
            else %transcending immunity
                DELTA(i,j)=1-exp(-contact(i,j)/delta); %delta^(contact(i,j));  %1-exp(-contact(i,j)/delta);  
            end
            
            if contact(i,j)==1
                NBOR(i,j)=epsilon; %mutation
            end
        end
    end
    
    %Count new infections
    newI=0;
    
    %Diff eqs
    dSdt = -(sum(immunity .* alpha .* S.*I/N)) + (gamma * sum(R));
    
    for i=1:n
        current_I=I(i);
        size(DELTA(i,:));
        size(I);
        size(R);
%         
%         mutate=.1;
%         dIdt(i)= alpha*S*I(i)/N - beta*I(i) + (sum(alpha * DELTA(i,:) .* R.*I(i) /N)); 
%         
%         nbor_to_mutate_to=randsample(1:length(contact(:,i)),1,true,contact(:,i));
%         dIdt(nbor_to_mutate_to)= dIdt(nbor_to_mutate_to) + mutate * alpha*S*I(i)/N;
        %%%%%dIdt(i)= alpha*S*I(i)/N - beta*I(i)+ (sum((NBOR(i,:).*I))) + (sum(alpha * DELTA(i,:) .* R.*I(i) /N)); 
        dIdt(i)= immunity(i) * alpha*S*I(i)/N - beta*I(i)+ (sum((NBOR(i,:).*I))) - (sum(NBOR(i,:).*I(i))) + (sum(alpha * DELTA(i,:) .* immunity(i) .* R.*I(i) /N)); 
        
        dRdt(i)= (beta * I(i)) - (gamma * R(i)) - ( alpha *R(i)* sum(DELTA(i,:).*I .* immunity/N ));
    end
    
    dSIR=[dSdt, dIdt, dRdt];
    dSIR=dSIR';
