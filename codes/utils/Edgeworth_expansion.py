from scipy.stats import norm
from scipy.special import hermitenorm
import numpy as np
from matplotlib import pyplot as plt


def phi(x, n=0):
    """计算\phi^n的值"""
    return ((-1)**n) * hermitenorm(n)(x) * norm.pdf(x, 0, 1)


def f_n(x, n, lambda_params):
    # 确保λ参数列表的长度正确

    if len(lambda_params) < 6:
        raise ValueError("lambda_params列表至少应包含6个元素")

    # 计算级数的各项
    term1 = phi(x)
    term2 = (1/n**0.5) * (1/6 * lambda_params[2] * phi(x, 3))
    term3 = (1/n) * (1/24 *
                     lambda_params[3] * phi(x, 4) + 1/72 * lambda_params[2]**2 * phi(x, 6))
    term4 = (1/n**1.5) * (1/120 * lambda_params[4] * phi(x, 5) + 1/144 * lambda_params[2]
                          * lambda_params[3] * phi(x, 7) + 1/1296 * lambda_params[2]**3 * phi(x, 9))
    term5 = (1/n**2) * (1/720 * lambda_params[5] * phi(x, 6) + (1/1152 * lambda_params[3]**2 + 1/720 * lambda_params[2] * lambda_params[4]) * phi(
        x, 8) + 1/1728 * lambda_params[2]**2 * lambda_params[3] * phi(x, 10) + 1/31104 * lambda_params[2]**4 * phi(x, 12))

    # 计算并返回函数的总值
    # return term1 - term2 + term3 - term4 + term5
    return [term1, term2, term3, term4, term5]

# 用来测试的λ参数列表
# lambda_params = [0, 0, 0.1, 0.1, 0.1, 0.1]  # 这些值仅仅是为了测试，您应该用实际的λ值来替换它们


def lambda_params_for_binomial(p, r):
    miu = p
    sigma = np.sqrt(p*(1-p))
    return p*((1-miu)/sigma)**r+(1-p)*((0-miu)/sigma)**r


# 调用函数并打印结果
p = 1e-6
lambda_params = [lambda_params_for_binomial(p, x) for x in range(7)]
plt.figure(figsize=(6, 4))

x_values = np.linspace(-15, 15, 100)  # 生成从-15到15之间的100个均匀分布的数值

terms = []
for x in x_values:
    term = f_n(x, 1e12, lambda_params)
    terms.append(term)

terms = np.array(terms).T  # 转置terms数组，使得每一列对应一个term
print(sum(abs(terms[0]))*30/99,sum(abs(terms[1]))*30/99,sum(abs(terms[2]))*30/99,sum(abs(terms[3]))*30/99,sum(abs(terms[4]))*30/99)

plt.semilogy(x_values, abs(terms[0]), label='term1')
plt.semilogy(x_values, abs(terms[1]), label='term2')
plt.semilogy(x_values, abs(terms[2]), label='term3')
plt.semilogy(x_values, abs(terms[3]), label='term4')
plt.semilogy(x_values, abs(terms[4]), label='term5')
plt.xlabel('x')
plt.ylabel('term')
plt.legend()
plt.show()
# plt.savefig('fig_edgeworth_expnsion.eps', dpi=600)
