from __future__ import division
import numpy as np
import lmoments


def dim_spi_n(values, scale_months, lower_limit=-3.09, upper_limit=3.09):
    '''
    :param values:
    :param scale_months:
    :param lower_limit:
    :param upper_limit:
    :return:
    '''

    #  replace the values array with sliding sums at the specified month scale
    values = get_sliding_sums(values, scale_months)

    # compute gamma parameters using the specified month scale
    gamma_values = gamma_parameters(values, scale_months)

    # replace the sums stored in the values array with fitted values
    probability = 0.0
    for month_index in range(scale_months - 1, len(values)):

        calendarMonth = month_index % 12
        if values[month_index] != np.nan:

            # compute the probability
            probability = lmoments.cdfgam(values[month_index], gamma_values[calendarMonth, :])

            # convert the probability to a fitted value
            values[month_index] = lmoments.quanor(probability, (0.0, 1.0))

    # return the fitted values clipped to the specified upper and lower limits
    return np.clip(values, lower_limit, upper_limit)


def gamma_parameters(summed_monthly_values, scale_months):
    '''
    :param monthly_values:
    :param scale_months:
    :return:
    '''

    # allocate the array of gamma parameters we'll return
    gamma_parameters = np.full((12, 2), np.nan, dtype=np.float64)

    # process each calendar month's values separately
    for i in range(12):

        # get the values for the calendar month
        calendar_month_sums = summed_monthly_values[i::12]

        # strip out all the NaN values
        calendar_month_sums = calendar_month_sums[np.logical_not(np.isnan(calendar_month_sums))]

        # get the non-zero values only (resulting array will still contain NaNs if present)
        nonzero_calendar_month_values = calendar_month_sums[np.nonzero(calendar_month_sums)]

        # determine how many zeros there were
        number_of_sums = calendar_month_sums.shape[0]
        number_of_nonzeros = nonzero_calendar_month_values.shape[0]
        number_of_zeros = number_of_sums - number_of_nonzeros

        # calculate the probability of zero, the first gamma parameter
        probability_of_zero = number_of_zeros / number_of_sums

        if probability_of_zero>0.1:
            gamma_parameters[i, :] = 0.0             
        else:
            # Fit gamma distribution
            try:
                LMU = lmoments.samlmu(nonzero_calendar_month_values)
                gamfit = lmoments.pelgam(LMU)
                gamma_parameters[i, :] = gamfit
            except:
                gamma_parameters[i, :] = 0.0

    return gamma_parameters


def get_sliding_sums(values, number_of_values_to_sum):
    '''
    Get the valid sliding summations using 1-D convolution. The initial (number_of_values_to_sum - 1) elements 
    of the result array will be padded with np.NaN values.

    :param values: the array of values over which we'll compute sliding sums
    :param number_of_values_to_sum: the number of values for which each sliding summation will encompass, for example if
            this value is 3 then the first two elements of the output array will contain the pad value and the third
            element of the output array will contain the sum of the first three elements, and so on
    :return: an array of sliding sums, equal in length to the input values array, left padded with NaN values
    '''
    # get the valid sliding summations with 1D convolution
    sliding_sums = np.convolve(values, np.ones(number_of_values_to_sum), mode='valid')

    # pad the first (n - 1) elements of the array with NaN values
    return np.hstack(([np.nan]*(number_of_values_to_sum - 1), sliding_sums))
