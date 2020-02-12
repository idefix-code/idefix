#include "idefix.hpp"
#include "input.hpp"
#include "gitversion.h"

Input::Input() {
    std::cout << "Creating Input\n";
    // Default constructor
    npoints[IDIR] = 64;
    npoints[JDIR] = 64;
    npoints[KDIR] = 64;

    xstart[IDIR] = 0.0;
    xstart[JDIR] = 0.0;
    xstart[KDIR] = 0.0;

    xend[IDIR] = 1.0;
    xend[JDIR] = 1.0;
    xend[KDIR] = 1.0;

    nghost[IDIR] = 2;
    nghost[JDIR] = 2;
    nghost[KDIR] = 2;

    // Number of integrator stages
    nstages=2;
    PrintLogo();
}

void Input::PrintLogo() {
    std::cout << "                                  .:HMMMMHn:.  ..:n.."<< std::endl;
    std::cout << "                                .H*'``     `'%HM'''''!x."<< std::endl;
    std::cout << "         :x                    x*`           .(MH:    `#h."<< std::endl;
    std::cout << "        x.`M                   M>        :nMMMMMMMh.     `n."<< std::endl;
    std::cout << "         *kXk..                XL  nnx:.XMMMMMMMMMMML   .. 4X."<< std::endl;
    std::cout << "          )MMMMMx              'M   `^?M*MMMMMMMMMMMM:HMMMHHMM."<< std::endl;
    std::cout << "          MMMMMMMX              ?k    'X ..'*MMMMMMM.#MMMMMMMMMx"<< std::endl;
    std::cout << "         XMMMMMMMX               4:    M:MhHxxHHHx`MMx`MMMMMMMMM>"<< std::endl;
    std::cout << "         XM!`   ?M                `x   4MM'`''``HHhMMX  'MMMMMMMM"<< std::endl;
    std::cout << "         4M      M                 `:   *>     `` .('MX   '*MMMM'"<< std::endl;
    std::cout << "          MX     `X.nnx..                        ..XMx`     'M*X"<< std::endl;
    std::cout << "           ?h.    ''```^'*!Hx.     :Mf     xHMh  M**MMM      4L`"<< std::endl;
    std::cout << "            `*Mx           `'*n.x. 4M>   :M` `` 'M    `       %"<< std::endl;
    std::cout << "             '%                ``*MHMX   X>      !"<< std::endl;
    std::cout << "            :!                    `#MM>  X>      `   :x"<< std::endl;
    std::cout << "           :M                        ?M  `X     .  ..'M"<< std::endl;
    std::cout << "           XX                       .!*X  `x   XM( MMx`h"<< std::endl;
    std::cout << "          'M>::                        `M: `+  MMX XMM `:"<< std::endl;
    std::cout << "          'M> M                         'X    'MMX ?MMk.Xx.."<< std::endl;
    std::cout << "          'M> ?L                     ...:!     MMX.H**'MMMM*h"<< std::endl;
    std::cout << "           M>  #L                  :!'`MM.    . X*`.xHMMMMMnMk."<< std::endl;
    std::cout << "           `!   #h.      :L           XM'*hxHMM*MhHMMMMMMMMMM'#h"<< std::endl;
    std::cout << "           +     XMh:    4!      x   :f   MM'   `*MMMMMMMMMM%  `X"<< std::endl;
    std::cout << "           M     Mf``tHhxHM      M>  4k xxX'      `#MMMMMMMf    `M .>"<< std::endl;
    std::cout << "          :f     M   `MMMMM:     M>   M!MMM:         '*MMf'     'MH*"<< std::endl;
    std::cout << "          !     Xf   'MMMMMX     `X   X>'h.`          :P*Mx.   .d*~.."<< std::endl;
    std::cout << "        :M      X     4MMMMM>     !   X~ `Mh.      .nHL..M#'%nnMhH!'`"<< std::endl;
    std::cout << "       XM      d>     'X`'**h     'h  M   ^'MMHH+*'`  ''''   `'**'"<< std::endl;
    std::cout << "    %nxM>      *x+x.:. XL.. `k     `::X"<< std::endl;
    std::cout << ":nMMHMMM:.  X>  Mn`*MMMMMHM: `:     ?MMn."<< std::endl;
    std::cout << "    `'**MML M>  'MMhMMMMMMMM  #      `M:^*x"<< std::endl;
    std::cout << "         ^*MMttnnMMMMMMMMMMMH>.        M:.4X"<< std::endl;
    std::cout << "                        `MMMM>X   (   .MMM:MM!   ."<< std::endl;
    std::cout << "                          `'''4x.dX  +^ `''MMMMHM?L.."<< std::endl;
    std::cout << "                                ``'           `'`'`'`"<< std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "           This is Idefix " << GITVERSION << std::endl;

}
