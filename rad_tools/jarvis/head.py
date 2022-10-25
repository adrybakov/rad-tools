from random import randint

from rad_tools.routines import GREEN, RESET


class Action:

    exits = {
        "exit",
        "exit()",
        "quit",
        "quit()",
        "kill",
        "shutdown"}
    tb2j = {
        "TB2J",
        "tb2j"
    }
    confusions = [
        "You`re confused?\nI`m fucking confused, bro!",
        "Sorry, I did not understand you, could you try again?",
        "Jarvis does not support this type of input, try better.",
    ]

    def __init__(self, name) -> None:
        self.name = name
        self.decide_what_to_do()

    def decide_what_to_do(self):
        if self.name in self.exits:
            return 0
        elif self.name in self.tb2j:
            print(f"{GREEN}Start to set up TB2J calculation{RESET}")
            tb2j_action = TB2JAction("Begin")
        elif self.name == "Initiate":
            pass
        else:
            self.confused()
        self.name = self.listen_to_user()

    def listen_to_user(self):
        message = input()
        self.name = message
        self.decide_what_to_do()

    def confused(self):
        i = randint(0, len(self.confusions) - 1)
        print(i)
        print(f"{GREEN}{self.confusions[i]}{RESET}")


class TB2JAction(Action):

    def __init__(self, name) -> None:
        super().__init__(name)
