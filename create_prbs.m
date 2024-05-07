function prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli)

%  CREATE_PRBS  is used for the generation of a PRBS signal
%
%        prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli)
%
%  "Entry parameters" are :
% 	  ValUinit  : Initial steady state
%     ValAmpli  : Magnitude
%     ValDecal  : Add-on DC component
%     ValLgReg  : Register length
%     ValDivi   : Frequency divider
%     Nsamp     : Number of samples
%     Tappli    : Application instant
%
%
%
%
% 	                ____  Valdecal + ValAmpli         __________      ____
%                  |    |                            |          |    |
%  Valdecal       -|----|--------                    |          |    |
%                  |    |____________________________|          |____|
%                  |
%                  |
%  ini ____________|
%                                                    |--------->|
%      |-Tappli -->|                        ValReg * ValDivi
%
%
%      |---------- Nsamp ------------------------------------------------->|
%
%
% 	"Exit parameter" is  :
%     prbs : prbs sequence created by PRBS algo


% initialization
k1 = ValLgReg - 1;
k2 = ValLgReg;

switch(ValLgReg)
    case 5
        k1 = 3;
    case 7
        k1 = 4;
    case 9
        k1 = 5;
    case 10
        k1 = 7;
    case 11
        k1 = 9;
end

for i=1:11
    sbpa(i) = 1;
end

% Run the PRBS algorithm

% Output set to init-value until the PRBS application istant

for i=1:Tappli
    prbs(i) = ValUinit;
end

% PRBS sequence generation
i=Tappli+1;
while (i<=Nsamp)
    uiu = -sbpa(k1)*sbpa(k2);
    if (ValLgReg == 8)
        uiu = -sbpa(2)*sbpa(3)*sbpa(5)*sbpa(8);
    end
    j=1;
    while (j<=ValDivi)
        prbs(i) = uiu * ValAmpli + ValDecal;
        i=i+1;
        j=j+1;
    end
    for j=ValLgReg:-1:2
        sbpa(j) = sbpa(j-1);
    end
    sbpa(1) = uiu;
end


