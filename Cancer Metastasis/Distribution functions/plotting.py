def plot_distribution_of_array(array , bins , the_title , x_lable , y_label):
    x=np.linspace(min(array), max(array), bins+1)
    y=[0]*(bins+1)
    for i in array:
        y[int((i-min(array))/(max(array)-min(array))*bins)]+=1
    plt.plot(x,np.array(y)/len(array) , marker = '^', linestyle = '--')
    plt.xlabel('{}'.format(x_lable))
    plt.ylabel('{}'.format(y_label))
    plt.title(the_title)
    plt.savefig('{}.png'.format(the_title))


def Plot_distribution_without_clf(array , bins , time ,Filename , type , N_cells):
    x_max = 1
    y_max=15
    x=np.linspace(min(array), max(array), bins+1)
    y=[0]*(bins+1)
    for i in array:
        y[int((i-min(array))/(max(array)-min(array))*bins)]+=1
    plt.scatter(x,np.array(y)/len(array) *100 )

    gauss = lambda x, sigma2_times_2_rev , mu: np.sqrt(sigma2_times_2_rev/np.pi) * np.exp(-1*sigma2_times_2_rev * (x-mu)**2)
    exponential = lambda x, a, b: a*np.exp(-b*x)
    percent = 0.06
    gauss_x = x[:int(percent*bins)]
    gauss_y = np.array(y)[:int(percent*bins)]/len(array) *100

    exp_x = x[int(percent*bins):]
    exp_y = np.array(y)[int(percent*bins):]/len(array) *100

    parameter_gauss, covarience_gauss = curve_fit(gauss, gauss_x, gauss_y , maxfev=5000)
    parameter_exp , covarience_exp = curve_fit(exponential, exp_x, exp_y , maxfev=5000)
    #parameter_gauss_all , covarience_gauss_all = curve_fit(gauss, x, np.array(y)/len(array) *100 , maxfev=5000)

    y_data_fit_gauss = gauss(gauss_x, parameter_gauss[0], parameter_gauss[1])
    y_data_fit_exp = exponential(exp_x, parameter_exp[0], parameter_exp[1])
    #y_data_fit_gauss_all = gauss(x, parameter_gauss_all[0], parameter_gauss_all[1])

    #plt.plot(gauss_x, y_data_fit_gauss, label='$\\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(\\Upsilon - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss[0]*2) ,parameter_gauss[1] ,1/(parameter_gauss[0]*2)) , color = 'g')
    plt.plot(gauss_x, y_data_fit_gauss, label='$ \\sigma^2 =  %5.3f , \\mu = %5.3f $'%( 1/(parameter_gauss[0]*2) , parameter_gauss[1]) , color = 'g')
    plt.plot(exp_x, y_data_fit_exp, label='$%5.3f e^{- %5.3f \\Upsilon} $' % (parameter_exp[0], parameter_exp[1]) , color = 'orange')
    #plt.plot(x, y_data_fit_gauss_all, label='$ fit: \\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(x - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss_all[0]*2) ,parameter_gauss_all[1] ,1/(parameter_gauss_all[0]*2)) , color = 'k')

    # if type == 'Soft':
    #     plt.scatter(x,np.array(y)/ N_cells * 100, label = type , color = 'blue')
    # elif type == 'Hard':
    #     plt.scatter(x,np.array(y)/ N_cells * 100, label = type , color = 'red')
    # else:
    #     plt.scatter(x,np.array(y)/ N_cells * 100, label = type , color = 'green')  

    plt.xlim(3.5 , 4.5)  
    plt.ylim(0.01,15)
    plt.yscale('log')   
    plt.xlabel('Cell Shape Parameter')
    plt.ylabel('N (%)')
    plt.legend()
    plt.title("Cell Shape parameter Distribution - time = {} ".format(time*10000))

