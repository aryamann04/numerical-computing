A = [0 1 0 0; 1 0 3 0; -0.5 0 -0.2 1; -0.5 -0.3 1 0]; 

P1 = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]; 
M1 = [1 0 0 0; 0 1 0 0; 0.5 0 1 0; 0.5 0 0 1]; 

P2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; 
M2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0.3 0 1]; 

P3 = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]; 
M3 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1.3/2.5 1]; 

% calculate ~M's  
m_1 = P3 * P2 * M1 * P2 * P3; 
m_2 = P3 * M2 * P3;
m_3 = M3;

m = m_3 * m_2 * m_1; % M~
P = P3 * P2 * P1; 
U = M3 * P3 * M2 * P2 * M1 * P1 * A;

L = inv(m)

PA = P * A;
LU = L * U;

norm(PA - LU) % check PA = LU

[L_matlab, U_matlab, P_matlab] = lu(A); 

norm(P_matlab - P)
norm(U_matlab - U)
norm(L_matlab - L)