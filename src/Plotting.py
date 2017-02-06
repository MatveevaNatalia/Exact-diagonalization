import matplotlib.pyplot as plt
import numpy as np

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

    ax.set_ylim(0, 1.05)

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

def single_level_plot(t1, t2, y1, y2, parameters):
    fig = plt.figure()
    ax = plt.subplot(111)

    ax.plot(t1, y1[:, 0], '-', color='black', label='1st level')
    ax.plot(t1, y1[:, 1], '-', color='red', label='2nd level')
    ax.plot(t1, y1[:, 2], '-', color='blue', label='3rd level')
    ax.plot(t1, y1[:, 3], '-', color='green', label='4th level')

    ax.plot(t2, y2[:, 0], '-', color='black')
    ax.plot(t2, y2[:, 1], '-', color='red')
    ax.plot(t2, y2[:, 2], '-', color='blue')
    ax.plot(t2, y2[:, 3], '-', color='green')


    #ax.ylim([0, :])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylim(0, 1.05)

    ax.set_xlabel('$t$', fontsize=18)
    ax.set_ylabel('$P_{{sl}}$', fontsize=18)



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

    pdfText = "nT_{0}_sl.pdf".format(parameters["num_T"])

    plt.savefig(pdfText, format="pdf", bbox_inches='tight')

    plt.show()

def current_plot(I_left_1, I_left_2, I_right_1, I_right_2, t1, t2, parameters, I_av_left, I_av_right):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(t1, I_left_1, '-', color='black', label='$I_L$')
    ax.plot(t1, I_right_1, '-', color='red', label='$I_R$')
    ax.plot(t2, I_left_2, '-', color='black')
    ax.plot(t2, I_right_2, '-', color='red')
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylim(bottom=0)

    ax.set_xlabel('$t$', fontsize=18)
    ax.set_ylabel('$I_L(t), I_R(t)$', fontsize=18)

    titleText = "$\gamma_R={0}\;, \gamma_L={1}\;, V_L={2}\;, V_R={3}$".format(
        parameters["gamma_R"],
        parameters["gamma_L"],
        parameters["V_L"],
        parameters["V_R"]
    )
    titleText += "\n $\mu_-={0}\;, \mu_+={1}\;, \Delta={2}\;,  V_{{sc}}={3}$".format(
        parameters["mu_minus"],
        parameters["mu_plus"],
        parameters["gap"],
        parameters["scat_ampl"]
    )

    plt.title(titleText, fontsize=16)

    labelText = "$T={0}\;, nT={1}\;,I_{{av}}^L={2:.2f}\;, I_{{av}}^R={3:.2f}\;$".format(
        parameters["T"],
        parameters["num_T"],
        I_av_left,
        I_av_right
    )

    plt.text(0.1, 0.9, labelText, transform=ax.transAxes, fontsize=16)

    pdfText = "current_{0}.pdf".format(parameters["num_T"])

    plt.savefig(pdfText, format="pdf", bbox_inches='tight')

    plt.show()

def current_average_plot(I_av_left, I_av_right, parameters):
    fig = plt.figure()
    ax = plt.subplot(111)

    num_T = parameters["num_T"]

    t = [i for i in range(num_T)]
    ax.plot(t, I_av_left, '-', color='black', label='$I_av_L$')
    ax.plot(t, I_av_right, '-', color='red', label='$I_av_R$')

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #ax.set_ylim(bottom=0)

    ax.set_xlabel('$t$', fontsize=18)
    ax.set_ylabel('$I_L(t), I_R(t)$', fontsize=18)

    titleText = "$\gamma_R={0}\;, \gamma_L={1}\;, V_L={2}\;, V_R={3}$".format(
        parameters["gamma_R"],
        parameters["gamma_L"],
        parameters["V_L"],
        parameters["V_R"]
    )
    titleText += "\n $\mu_-={0}\;, \mu_+={1}\;, \Delta={2}\;,  V_{{sc}}={3}$".format(
        parameters["mu_minus"],
        parameters["mu_plus"],
        parameters["gap"],
        parameters["scat_ampl"]
    )

    plt.title(titleText, fontsize=16)

    labelText = "$T={0}\;, nT={1}\;$".format(
        parameters["T"],
        parameters["num_T"],
        I_av_left,
        I_av_right
    )

    plt.text(0.1, 0.9, labelText, transform=ax.transAxes, fontsize=16)

    pdfText = "current_average.pdf"

    plt.savefig(pdfText, format="pdf", bbox_inches='tight')

    plt.show()

def current_average_plot_file(file_name, parameters):
    fig = plt.figure()
    ax = plt.subplot(111)

    data = np.loadtxt(file_name)

    ax.plot(data[:, 0], data[:, 1], 'o', ms = 6, color='black',
            markeredgecolor='black', label='$I_{{av}}^L$')
    ax.plot(data[:, 0], data[:, 2], 's', ms=6, color='red',
            markeredgecolor='red', label='$I_{{av}}^R$')

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints=1)

    #ax.set_ylim(bottom=0)

    ax.set_xlabel('$t/T$', fontsize=18)
    ax.set_ylabel('$I_{{av}}^L, I_{{av}}^R$', fontsize=18)

    titleText = "$\gamma_R={0}\;, \gamma_L={1}\;, V_L={2}\;, V_R={3}$".format(
        parameters["gamma_R"],
        parameters["gamma_L"],
        parameters["V_L"],
        parameters["V_R"]
    )
    titleText += "\n $\mu_-={0}\;, \mu_+={1}\;, \Delta={2}\;,  V_{{sc}}={3}\;, N_{{points}}={4}$".format(
        parameters["mu_minus"],
        parameters["mu_plus"],
        parameters["gap"],
        parameters["scat_ampl"],
        parameters["num_points"]
    )

    plt.title(titleText, fontsize=16)

    pdfText = "current_average.pdf"

    plt.savefig(pdfText, format="pdf", bbox_inches='tight')

    plt.show()

