% pre burner parameters
m_z = 685;
R = 3.6;
At = 0.0426;
Pf1 = 0.4;
Po2 = 0.4;
Ppbf = 61.1;

%intial Parameters
Rpb = [Rpbf ; Rpbo]; %3.6
k = 0;



%MwO = ; 
Hf = -89233.000;
MwF = 16.4;
Fu = Hf/MwF;

%w0 = ;
%Rc2c = ;

f = @(Rpbo) ((Rc2*Rpbo)/(1 + Rpbo) - R*w0);
df =  @(x) -0.5 * x.^(-3/2) - ((2.51/Re) ./ (r/3.7 + 2.51 ./ (Re*sqrt(x)))).*x.^(-1.5);

