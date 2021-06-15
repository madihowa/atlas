import matplotlib.pyplot as plt


def plt_result(inputs, predictions, target):


        fig = plt.figure()
        ax=fig.add_axes([0,0,1,1])
        ax.scatter(predictions, target, color='b')
        ax.plot(predictions, predictions, color='r')
        ax.set_xlabel('Grades Range')
        ax.set_ylabel('Grades Scored')
        ax.set_title('scatter plot')
        plt.xlabel('True Energy')
        plt.xlabel('Calibrated energy')
        #plt.show()
