import matplotlib.pyplot as plt

#def Labels(size_Np1, size_N): # To make it friend of tunneling ???
#    q = [0] * size_Np1
#    p = [0] * size_N

#    for i in range(size_N):
#        p[i] = str(to_bitfield(H_data_N.basisFock[i][1], num_level))

#        for i in range(H_data_Np1.size):
#        q[i] = str(to_bitfield(H_data_Np1.basisFock[i][1], num_level))


def populations_plot(t1, t2, y1, y2, p, q, parameters):
    fig = plt.figure()
    ax = plt.subplot(111)

    ax.plot(t1, y1[:, 0], '--', color='black', label=p[0])
    ax.plot(t1, y1[:, 1], '--', color='red', label=p[1])
    ax.plot(t1, y1[:, 2], '--', color='blue', label=p[2])
    ax.plot(t1, y1[:, 3], '--', color='green', label=p[3])
    ax.plot(t1, y1[:, 4], '-', color='black', label=q[0])
    ax.plot(t1, y1[:, 5], '-', color='red', label=q[1])
    ax.plot(t1, y1[:, 6], '-', color='blue', label=q[2])
    ax.plot(t1, y1[:, 7], '-', color='green', label=q[3])
    ax.plot(t1, y1[:, 8], '-', color='black', label=q[4])
    ax.plot(t1, y1[:, 9], '-', color='red', label=q[5])
    ax.plot(t2, y2[:, 0], '--', color='black')
    ax.plot(t2, y2[:, 1], '--', color='red')
    ax.plot(t2, y2[:, 2], '--', color='blue')
    ax.plot(t2, y2[:, 3], '--', color='green')
    ax.plot(t2, y2[:, 4], '-', color='black')
    ax.plot(t2, y2[:, 5], '-', color='red')
    ax.plot(t2, y2[:, 6], '-', color='blue')
    ax.plot(t2, y2[:, 7], '-', color='green')
    ax.plot(t2, y2[:, 8], '-', color='black')
    ax.plot(t2, y2[:, 9], '-', color='red')

    #ax.ylim([0, :])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylim(bottom=0)

    ax.set_xlabel('$t$', fontsize=18)
    ax.set_ylabel('$p(t), q(t)$', fontsize=18)



    titleText  = "$\gamma_R={0}\;, \gamma_L={1}\;, V_L={2}\;, V_R={3}$".format(
        parameters["gamma_R"],
        parameters["gamma_L"],
        parameters["V_L"],
        parameters["V_R"]
    )
    titleText += "\n $\mu_-={0}\;, \mu_+={1}\;, \Delta={2}\;,  V_{{sc}}={3}$ ".format(
        parameters["mu_minus"],
        parameters["mu_plus"],
        parameters["gap"],
        parameters["scat_ampl"]
    )

    plt.title(titleText, fontsize=16)

    labelText = "$T={0}\;, nT={1}$".format(
        parameters["T"],
        parameters["num_T"]
    )

    plt.text(0.1, 0.9, labelText, transform = ax.transAxes, fontsize=16)

    pdfText = "nT_{0}.pdf".format(parameters["num_T"])

    plt.savefig(pdfText, format="pdf", bbox_inches='tight')

    plt.show()

