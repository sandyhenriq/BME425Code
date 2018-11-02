function [ah, bh, aj, bj,am, bm, ad, bd, af, bf, ...
    ax, bx, ak1, bk1, xigate ] = LR1( Vm, EK1 )

am = 0.32 .* (Vm + 47.13) ./ (1 - exp(-.1 .* (Vm + 47.13)));
bm = 0.08 .* exp(-Vm./11);
ad = 0.095 .* exp(-.01 .* (Vm - 5)) ./ (1 + exp(-.072 .* (Vm - 5)));
bd = 0.07 .* exp(-.017 .* (Vm + 44)) ./ (1 + exp(.05 .* (Vm + 44)));
af = 0.012 .* exp(-.008 .* (Vm + 28)) ./ (1 + exp(.15 .* (Vm + 28)));
bf = 0.0065 .* exp(-.02 .* (Vm + 30)) ./ (1 + exp(-.2 .* (Vm + 30)));
ax = 0.0005 .* exp(.083 .* (Vm + 50)) ./ (1 + exp(.057 .* (Vm + 50)));
bx = .0013 .* exp(-.06 .* (Vm + 20)) ./ (1 + exp(-.04 .* (Vm + 20)));
ak1 = 1.02 ./ (1 + exp(.2385 .* (Vm - EK1 - 59.215)));
bk1 = (.49124 .* exp(.08032 .* (Vm - EK1 + 5.476)) + exp(.06175 .* ...
    (Vm - EK1 - 594.31))) ./ (1 + exp(-.5143 .* (Vm - EK1 + 4.753)));

for i=1:length(Vm)
if Vm(i) >= -40
    ah(i) = 0;
    aj(i) = 0;
    bh(i) = 1 ./ (0.13 .* (1 + exp((Vm(i) + 10.66) ./ -11.1)));
    bj(i) = 0.3 .* exp(-.0000002535) ./ (1 + exp(-.1 .* (Vm(i) + 32)));
else
    ah(i) = .135 .* exp((80 + Vm(i)) ./ -6.8);
    bh(i) = 3.56 .* exp(0.079 .* Vm(i)) + .000031 .* exp(0.35 .* Vm(i));
    aj(i) = (-.000012714 .* exp(.2444 .* Vm(i)) - .00003474 .* exp(-.04391 .* Vm(i))) ...
        .* (Vm(i) + 37.78) ./ (1 + exp(0.311 .* (Vm(i) + 79.23)));
    bj(i) = 0.1212 .* exp(-.01052 .* Vm(i)) ./ (1 + exp(-0.1378 .* (Vm(i) + 40.14)));
end
    if Vm(i) > -100
        xigate(i) = 2.837 .* (exp(0.04 .* (Vm(i) + 77)) - 1) ./ ((Vm(i) + 77) .* ... 
            exp(0.04 .* (Vm(i) + 35)));
    else
        xigate(i) = 1;
    end
end
end
