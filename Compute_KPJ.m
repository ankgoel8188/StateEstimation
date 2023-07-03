function [K,P,J] = Compute_KPJ(System, P, Type )


A = System.A;
% B = System.B;
C = System.C;

Q = System.Q;
R = System.R;
n = size(A,1);
if isfield(System,'Gamma')
    Gamma = System.Gamma;
else
    Gamma = eye(n);
end


pi_k    = Gamma * (Gamma' * Gamma)^(-1) * Gamma' ;
pi_k_p  = eye(n) - pi_k;
Impi = pi_k_p;
% Impi = System.Impi;


if sum(sum(abs(P)>1e10)) && 0 
    K = zeros(size(Gamma,2), size(C,1));
    P = P;
    J = trace(P);
    
else    
    warning('');
    if strcmp(Type, 'OP')
        K = inv(Gamma' * Gamma) * Gamma' * A * P * C' * inv(C * P * C' + R);
        
        Pm1 = P;
        
        P = A*P*A' + Q*eye(n) + ...
            (Gamma * K) * (C * P * C' + R) *(Gamma * K)' - ...
            A * P * C' * K' * Gamma' - ...
            (A * P * C' * K' * Gamma')';
        
%         MT = (A*Pm1*C')*inv(C*Pm1*C'+R)*(A*Pm1*C')';
%         P1 = A*Pm1*A' + Q*eye(n) ...
%             - MT ...
%             + Impi * MT * Impi';
% %         if norm(P-P1)
% %             disp(norm(P-P1) )
% %         end
%         
%         if norm(P1-P1')
%             disp(norm(P1-P1'))
%         end
        
        J = trace(P);
        
    end
    
    if strcmp(Type, 'RC')
        K = System.K;
        P = A*P*A' + Q*eye(n) + ...
            (Gamma * K) * (C * P * C' + R) *(Gamma * K)' - ...
            A * P * C' * K' * Gamma' - ...
            (A * P * C' * K' * Gamma')';
        J = trace(P);
    end
    
    if strcmp(Type, 'RCF')
        K = System.K;
        %         P = A*P*A' + Q*eye(n) + ...
        %             (Gamma * K) * (C * P * C' + R) *(Gamma * K)' - ...
        %             A * P * C' * K' * Gamma' - ...
        %             (A * P * C' * K' * Gamma')';
        Pbar = A*P*A'+Q*eye(n);
        Rbar = C*Pbar*C' + R;
        P = A*P*A' + Q*eye(n) + ...
            (Gamma * K) * Rbar *(Gamma * K)' - ...
            Pbar * C' * K' * Gamma' - ...
            (Pbar * C' * K' * Gamma')';
        J = trace(P);
    end
    
    if strcmp(Type, 'NoGain')
        K = System.K*0;
        P = A*P*A' + Q*eye(n);
        J = trace(P);
    end
    
    
    if strcmp(Type, 'OF')
        Pbar = A*P*A'+Q*eye(n);
        Rbar = C*Pbar*C' + R;
        K = inv(Gamma' * Gamma) * Gamma' * Pbar * C' * inv(Rbar);
        P = A*P*A' + Q*eye(n) + ...
            (Gamma * K) * Rbar *(Gamma * K)' - ...
            Pbar * C' * K' * Gamma' - ...
            (Pbar * C' * K' * Gamma')';
        
        J = trace(P);
        
    end
    
     if strcmp(Type, 'Filt_K')
        Pbar = A*P*A'+Q;
        Rbar = C*Pbar*C' + R;
        K = System.K;;
        P = A*P*A' + Q*eye(n) + ...
            (Gamma * K) * Rbar *(Gamma * K)' - ...
            Pbar * C' * K' * Gamma' - ...
            (Pbar * C' * K' * Gamma')';
        
        J = trace(P);
    end
    
%     [warnMsg, warnId] = lastwarn;
%     if ~isempty(warnMsg)
%         keyboard
%     end
    
    
    
    
    if strcmp(Type, '1SP')
        K = inv(Gamma' * Gamma) * Gamma' * A * P * C' * inv(C * P * C' + R);
        P = A*P*A' + Q*eye(n) + ...
            Gamma * K * (C * P * C' + R) *(Gamma * K)' - ...
            A * P * C' * K' * Gamma' - ...
            (A * P * C' * K' * Gamma')';
        J = trace(P);
    end
    
    if strcmp(Type, '2SP')
        K = inv(Gamma' * A' * A * Gamma) * Gamma' * A' * A * P * C' * inv(C * P * C' + R);
        P = A*P*A' + Q*eye(n) + ...
            A * Gamma * K * (C * P * C' + R) *(A * Gamma * K)' - ...
            A * P * C' * K' * Gamma' * A' - ...
            A * (P * C' * K' * Gamma')' * A';
        J = trace(P);
    end
    
    if strcmp(Type, '2SF')
        K = inv(Gamma' * Gamma) * Gamma' * (A * P * A' + Q) * C' * inv(C * P * C' + R);
        P = (eye(n) - K * C) * P;
        J = trace(P);
    end
    
    
    if strcmp(Type, '1SP_P_only')
        K = System.K;
        P = A*P*A' + Q*eye(n) + ...
            Gamma * K * (C * P * C' + R) *(Gamma * K)' - ...
            A * P * C' * K' * Gamma' - ...
            (A * P * C' * K' * Gamma')';
        J = trace(P);
    end
    
    if strcmp(Type, '2SP_P_only')
        K = System.K;
        P = A*P*A' + Q*eye(n) + ...
            A * Gamma * K * (C * P * C' + R) *(A * Gamma * K)' - ...
            A * P * C' * K' * Gamma' * A' - ...
            A * (P * C' * K' * Gamma')' * A';
        J = trace(P);
    end
    
end

end

