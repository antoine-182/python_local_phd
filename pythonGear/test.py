# -*- coding: utf-8 -*-

# import argparse
# # https://docs.python.org/3/library/argparse.html
# parser = argparse.ArgumentParser(description='Add some integers.')
#
# # argument positionnel
# parser.add_argument('integers', metavar='N', type=int,
#                     nargs='+', # iterable
#                     help='integer list')
# # argument optionnel
# # action='store_const' = si rien n'est précisé, active tout de même ?
# # dest : nom qu'on lui donne, autrement il s'appelerait sum
# parser.add_argument('-s','--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
#
#
# args = parser.parse_args()
# print(args.accumulate(args.integers))

""" https://docs.python.org/fr/3/howto/argparse.html
    Premier exemple parlant
"""
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("square", type=int,
#                     help="display a square of a given number")
# parser.add_argument("-v", "--verbose", action="store_true",
#                     help="increase output verbosity")
# args = parser.parse_args()
# answer = args.square**2
# if args.verbose:
#     print("the square of {} equals {}".format(args.square, answer))
# else:
#     print(answer)

""" https://docs.python.org/fr/3/howto/argparse.html
    Deuxième exemple, plus pratique
"""
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("x", type=int, help="the base")
# parser.add_argument("y", type=int, help="the exponent")
# parser.add_argument("-v", "--verbosity", action="count", default=0)
# args = parser.parse_args()
# answer = args.x**args.y
# if args.verbosity >= 2:
#     print("Running '{}'".format(__file__))
# if args.verbosity >= 1:
#     print("{}^{} == ".format(args.x, args.y), end="")
# print(answer)

""" https://docs.python.org/fr/3/howto/argparse.html
    Troisième exemple conflit
"""
# import argparse
#
# parser = argparse.ArgumentParser(description="calculate X to the power of Y")
# group = parser.add_mutually_exclusive_group()
# group.add_argument("-v", "--verbose", action="store_true")
# group.add_argument("-q", "--quiet", action="store_true")
# parser.add_argument("x", type=int, help="the base")
# parser.add_argument("y", type=int, help="the exponent")
# args = parser.parse_args()
# answer = args.x**args.y
#
# if args.quiet:
#     print(answer)
# elif args.verbose:
#     print("{} to the power {} equals {}".format(args.x, args.y, answer))
# else:
#     print("{}^{} == {}".format(args.x, args.y, answer))

""" Mon exemple
"""
