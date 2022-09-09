function a=propensity_for_RB(X,k,phi,reactionI)
% The propensity function of Linear_metabolite_pathway for RB simulations

a=zeros(1,size(X,2));

a=a+(k(1)).*(1==reactionI);
a=a+(k(2).*X(1,:)).*(2==reactionI);
a=a+(k(3).*X(1,:)).*(3==reactionI);
a=a+(k(4)).*(4==reactionI);
a=a+(k(5).*X(2,:)).*(5==reactionI);
a=a+(k(6).*X(2,:)).*(6==reactionI);
a=a+(k(7)).*(7==reactionI);
a=a+(k(8).*X(3,:)).*(8==reactionI);
a=a+(k(9).*X(3,:)).*(9==reactionI);
a=a+(k(10)).*(10==reactionI);
a=a+(k(11).*X(4,:)).*(11==reactionI);
a=a+(k(12).*X(4,:)).*(12==reactionI);
a=a+(k(13)).*(13==reactionI);
a=a+(k(14).*X(5,:)).*(14==reactionI);
a=a+(k(15).*X(5,:)).*(15==reactionI);
a=a+(k(16)).*(16==reactionI);
a=a+(k(17).*X(6,:)).*(17==reactionI);
a=a+(k(18).*X(6,:)).*(18==reactionI);
a=a+(k(19)).*(19==reactionI);
a=a+(k(20).*X(7,:)).*(20==reactionI);
a=a+(k(21).*X(7,:)).*(21==reactionI);
a=a+(k(22)).*(22==reactionI);
a=a+(k(23).*X(8,:)).*(23==reactionI);
a=a+(k(24).*X(8,:)).*(24==reactionI);
a=a+(k(25)).*(25==reactionI);
a=a+(k(26).*X(9,:)).*(26==reactionI);
