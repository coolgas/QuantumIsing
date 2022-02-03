import numpy as np
import matplotlib.pyplot as plt
import os


class CountCalls:
    def __init__(self, func):
        self._count = 0
        self._func = func

    def __call__(self, *args, **kwargs):
        self._count += 1
        return self._func(*args, **kwargs)

    @property
    def call_count(self):
        return self._count



@CountCalls
def window_seeker(last_day_noise_diff, stock_picker_name, lst_x, lst_y, n_seq=0, steps=0, pos_cut=0, neg_cut=0, mean_first=True, mean_cut=0, std_cut=0):
    '''
    Inputs:
        n_seq: the number of how many segmentations you want to seperate the domain into.
        pos_cut: the cut of the number of postive y's in the domain.
        neg_cut: the cut of the number of negative y's in the domain.
        step: how many steps you want to scan through x.
        mean_first: if set to true, the program will extract those domains with the highest means.
        mean_cut: the cut of mean value below which don't need to be considered.
        std_cut: the cut of standard deviation below which don't need to be considered.
    '''
    domain_len = (lst_x[-1] - lst_x[0]) / n_seq
    steps_len = (lst_x[-1] - lst_x[0]) / steps
    
    mean_lst = [] # this list would store the mean value of each domain
    std_lst = [] # this list would store the std value of each domain
    pos_neg_lst = [] # this list would store the number of postive and negative y's as a pair in each domain

    xy_merged = list(zip(lst_x, lst_y))
    
    i = 0
        
    while i <= steps - 1:
        x0 = lst_x[0] + i * steps_len # the beginning point of the domain
        x1 = lst_x[0] + domain_len + i * steps_len # the ending point of the domain
        
        domain_y_lst = [] # this list will store the y in the domain for std calculations

        y_total = 0
        pos_y_num = 0 # the number of positive y's
        neg_y_num = 0 # the number of negative y's
        
        for xy_pair in xy_merged:
            x = xy_pair[0]
            y = xy_pair[1]
            if x0 <= x <= x1:
                domain_y_lst.append(y)
                y_total += y
                if y > 0:
                    pos_y_num += 1
                else:
                    neg_y_num += 1
        
        pos_neg_lst.append((pos_y_num, neg_y_num))

        # the only condition we would consider
        if pos_y_num >= pos_cut and neg_y_num <= neg_cut:
            mean_lst.append(y_total/y_num)
            std_lst.append(np.std(domain_y_lst))
        else:
            mean_lst.append(0)
            std_lst.append(10**20) # this plays the role as place holder 
        
        i += 1
    
    if mean_first:
        # extract the maximum mean from the mean_lst
        max_mean = max(mean_lst)
        max_index = mean_lst.index(max_mean)
        if max_mean >= mean_cut:
            x0 = lst_x[0] + max_index * steps_len
            x1 = lst_x[0] + domain_len + max_index * steps_len

            print("GOLDEN STICK FOUND!!!")
            print("we are in the mode: \'mean_first\'")
            print("label: {}".format(window_seeker.call_count))
            print("n_seg: {}".format(n_seg))
            print("domain_mean: {}".format(max_mean))
            print("domain_std: {}".format(std_lst[max_ind]))
            print("n_pos: {}".format(pos_neg_lst[max_index][0]))
            print("n_neg: {}".format(pos_neg_lst[max_index][1]))
            print("domain: [{}, {}]".format(x0, x1))

            if not os.path.exists('Windows'):
                os.makedirs('Windows')

            fig = plt.figure(figsize=(10, 5))
            fig.patch.set_facecolor('xkcd:silver')
            plt.scatter(lst_x, lst_y)
            plt.title(stock_picker_name, label='label: {}, domain: [{}, {}]'.format(window_seeker.call_count, x0, x1))
            plt.axhline(y=0.0, color='black', linestyle='--')
            plt.axvline(x=x0, color='red', linestyle='--')
            plt.axvline(x=x1, color='red', linestyle='--')
            plt.plot(last_day_noise_diff, 0 , 'x', color='red')
            plt.close(fig)
            fig.savefig('./Windows/' + stock_picker_name + '.png')
            return [x0, x1]

        else:
            print("NO GOLDEN STICK FOR THIS GRAPH!!!")
            print("=================================================")

    else:
        # extract the minimum std from the std_lst
        min_std = min(std_lst)
        min_index = std_lst.index(min_std)
        if min_std <= std_cut:
            x0 = lst_x[0] + min_index * steps_len
            x1 = lst_x[0] + domain_len + min_index * steps_len

            print("GOLDEN STICK FOUND!!!")
            print("we are in the mode: \'std_first\'")
            print("label: {}".format(window_seeker.call_count))
            print("n_seg: {}".format(n_seg))
            print("domain_mean: {}".format(mean_lst[min_index]))
            print("domain_std: {}".format(min_std))
            print("n_pos: {}".format(pos_neg_lst[min_index][0]))
            print("n_neg: {}".format(pos_neg_lst[min_index][1]))
            print("domain: [{}, {}]".format(x0, x1))

            if not os.path.exists('Windows'):
                os.makedirs('Windows')

            fig = plt.figure(figsize=(10, 5))
            fig.patch.set_facecolor('xkcd:silver')
            plt.scatter(lst_x, lst_y)
            plt.title(stock_picker_name, label='label: {}, domain: [{}, {}]'.format(window_seeker.call_count, x0, x1))
            plt.axhline(y=0.0, color='black', linestyle='--')
            plt.axvline(x=x0, color='red', linestyle='--')
            plt.axvline(x=x1, color='red', linestyle='--')
            plt.plot(last_day_noise_diff, 0 , 'x', color='red')
            plt.close(fig)
            fig.savefig('./Windows/' + stock_picker_name + '.png')
            return [x0, x1]

        else:
            print("NO GOLDEN STICK FOR THIS GRAPH!!!")
            print("=================================================")

    return


        


