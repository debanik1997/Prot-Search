import React from "react";
import { BrowserRouter as Router, Switch, Route } from "react-router-dom";
import SearchPage from "./SearchPage";
import Results from "./Results";
function App() {
  return (
    <Router>
      <Switch>
        <Route
          exact
          path="/results"
          component={(props) => <Results {...props} />}
        />
        <Route path="/" component={SearchPage} />
        <Route path="*" component={SearchPage} />
      </Switch>
    </Router>
  );
}

export default App;
